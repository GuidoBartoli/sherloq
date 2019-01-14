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
#include "alglibinternal.h"

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


}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS IMPLEMENTATION OF COMPUTATIONAL CORE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{




static void tsort_tagsortfastirec(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Integer */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);
static void tsort_tagsortfastrrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);
static void tsort_tagsortfastrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* bufa,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);




















static void hsschur_internalauxschur(ae_bool wantt,
     ae_bool wantz,
     ae_int_t n,
     ae_int_t ilo,
     ae_int_t ihi,
     /* Real    */ ae_matrix* h,
     /* Real    */ ae_vector* wr,
     /* Real    */ ae_vector* wi,
     ae_int_t iloz,
     ae_int_t ihiz,
     /* Real    */ ae_matrix* z,
     /* Real    */ ae_vector* work,
     /* Real    */ ae_vector* workv3,
     /* Real    */ ae_vector* workc1,
     /* Real    */ ae_vector* works1,
     ae_int_t* info,
     ae_state *_state);
static void hsschur_aux2x2schur(double* a,
     double* b,
     double* c,
     double* d,
     double* rt1r,
     double* rt1i,
     double* rt2r,
     double* rt2i,
     double* cs,
     double* sn,
     ae_state *_state);
static double hsschur_extschursign(double a, double b, ae_state *_state);
static ae_int_t hsschur_extschursigntoone(double b, ae_state *_state);




static ae_bool safesolve_cbasicsolveandupdate(ae_complex alpha,
     ae_complex beta,
     double lnmax,
     double bnorm,
     double maxgrowth,
     double* xnorm,
     ae_complex* x,
     ae_state *_state);


static ae_bool hpccores_hpcpreparechunkedgradientx(/* Real    */ ae_vector* weights,
     ae_int_t wcount,
     /* Real    */ ae_vector* hpcbuf,
     ae_state *_state);
static ae_bool hpccores_hpcfinalizechunkedgradientx(/* Real    */ ae_vector* buf,
     ae_int_t wcount,
     /* Real    */ ae_vector* grad,
     ae_state *_state);


static void xblas_xsum(/* Real    */ ae_vector* w,
     double mx,
     ae_int_t n,
     double* r,
     double* rerr,
     ae_state *_state);
static double xblas_xfastpow(double r, ae_int_t n, ae_state *_state);


static double linmin_ftol = 0.001;
static double linmin_xtol = 100*ae_machineepsilon;
static ae_int_t linmin_maxfev = 20;
static double linmin_stpmin = 1.0E-50;
static double linmin_defstpmax = 1.0E+50;
static double linmin_armijofactor = 1.3;
static void linmin_mcstep(double* stx,
     double* fx,
     double* dx,
     double* sty,
     double* fy,
     double* dy,
     double* stp,
     double fp,
     double dp,
     ae_bool* brackt,
     double stmin,
     double stmax,
     ae_int_t* info,
     ae_state *_state);


static ae_bool ntheory_isprime(ae_int_t n, ae_state *_state);
static ae_int_t ntheory_modmul(ae_int_t a,
     ae_int_t b,
     ae_int_t n,
     ae_state *_state);
static ae_int_t ntheory_modexp(ae_int_t a,
     ae_int_t b,
     ae_int_t n,
     ae_state *_state);


static ae_int_t ftbase_coltype = 0;
static ae_int_t ftbase_coloperandscnt = 1;
static ae_int_t ftbase_coloperandsize = 2;
static ae_int_t ftbase_colmicrovectorsize = 3;
static ae_int_t ftbase_colparam0 = 4;
static ae_int_t ftbase_colparam1 = 5;
static ae_int_t ftbase_colparam2 = 6;
static ae_int_t ftbase_colparam3 = 7;
static ae_int_t ftbase_colscnt = 8;
static ae_int_t ftbase_opend = 0;
static ae_int_t ftbase_opcomplexreffft = 1;
static ae_int_t ftbase_opbluesteinsfft = 2;
static ae_int_t ftbase_opcomplexcodeletfft = 3;
static ae_int_t ftbase_opcomplexcodelettwfft = 4;
static ae_int_t ftbase_opradersfft = 5;
static ae_int_t ftbase_opcomplextranspose = -1;
static ae_int_t ftbase_opcomplexfftfactors = -2;
static ae_int_t ftbase_opstart = -3;
static ae_int_t ftbase_opjmp = -4;
static ae_int_t ftbase_opparallelcall = -5;
static ae_int_t ftbase_maxradix = 6;
static ae_int_t ftbase_updatetw = 16;
static ae_int_t ftbase_recursivethreshold = 1024;
static ae_int_t ftbase_raderthreshold = 19;
static ae_int_t ftbase_ftbasecodeletrecommended = 5;
static double ftbase_ftbaseinefficiencyfactor = 1.3;
static ae_int_t ftbase_ftbasemaxsmoothfactor = 5;
static void ftbase_ftdeterminespacerequirements(ae_int_t n,
     ae_int_t* precrsize,
     ae_int_t* precisize,
     ae_state *_state);
static void ftbase_ftcomplexfftplanrec(ae_int_t n,
     ae_int_t k,
     ae_bool childplan,
     ae_bool topmostplan,
     ae_int_t* rowptr,
     ae_int_t* bluesteinsize,
     ae_int_t* precrptr,
     ae_int_t* preciptr,
     fasttransformplan* plan,
     ae_state *_state);
static void ftbase_ftpushentry(fasttransformplan* plan,
     ae_int_t* rowptr,
     ae_int_t etype,
     ae_int_t eopcnt,
     ae_int_t eopsize,
     ae_int_t emcvsize,
     ae_int_t eparam0,
     ae_state *_state);
static void ftbase_ftpushentry2(fasttransformplan* plan,
     ae_int_t* rowptr,
     ae_int_t etype,
     ae_int_t eopcnt,
     ae_int_t eopsize,
     ae_int_t emcvsize,
     ae_int_t eparam0,
     ae_int_t eparam1,
     ae_state *_state);
static void ftbase_ftpushentry4(fasttransformplan* plan,
     ae_int_t* rowptr,
     ae_int_t etype,
     ae_int_t eopcnt,
     ae_int_t eopsize,
     ae_int_t emcvsize,
     ae_int_t eparam0,
     ae_int_t eparam1,
     ae_int_t eparam2,
     ae_int_t eparam3,
     ae_state *_state);
static void ftbase_ftapplysubplan(fasttransformplan* plan,
     ae_int_t subplan,
     /* Real    */ ae_vector* a,
     ae_int_t abase,
     ae_int_t aoffset,
     /* Real    */ ae_vector* buf,
     ae_int_t repcnt,
     ae_state *_state);
static void ftbase_ftapplycomplexreffft(/* Real    */ ae_vector* a,
     ae_int_t offs,
     ae_int_t operandscnt,
     ae_int_t operandsize,
     ae_int_t microvectorsize,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
static void ftbase_ftapplycomplexcodeletfft(/* Real    */ ae_vector* a,
     ae_int_t offs,
     ae_int_t operandscnt,
     ae_int_t operandsize,
     ae_int_t microvectorsize,
     ae_state *_state);
static void ftbase_ftapplycomplexcodelettwfft(/* Real    */ ae_vector* a,
     ae_int_t offs,
     ae_int_t operandscnt,
     ae_int_t operandsize,
     ae_int_t microvectorsize,
     ae_state *_state);
static void ftbase_ftprecomputebluesteinsfft(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* precr,
     ae_int_t offs,
     ae_state *_state);
static void ftbase_ftbluesteinsfft(fasttransformplan* plan,
     /* Real    */ ae_vector* a,
     ae_int_t abase,
     ae_int_t aoffset,
     ae_int_t operandscnt,
     ae_int_t n,
     ae_int_t m,
     ae_int_t precoffs,
     ae_int_t subplan,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     /* Real    */ ae_vector* bufc,
     /* Real    */ ae_vector* bufd,
     ae_state *_state);
static void ftbase_ftprecomputeradersfft(ae_int_t n,
     ae_int_t rq,
     ae_int_t riq,
     /* Real    */ ae_vector* precr,
     ae_int_t offs,
     ae_state *_state);
static void ftbase_ftradersfft(fasttransformplan* plan,
     /* Real    */ ae_vector* a,
     ae_int_t abase,
     ae_int_t aoffset,
     ae_int_t operandscnt,
     ae_int_t n,
     ae_int_t subplan,
     ae_int_t rq,
     ae_int_t riq,
     ae_int_t precoffs,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
static void ftbase_ftfactorize(ae_int_t n,
     ae_bool isroot,
     ae_int_t* n1,
     ae_int_t* n2,
     ae_state *_state);
static ae_int_t ftbase_ftoptimisticestimate(ae_int_t n, ae_state *_state);
static void ftbase_ffttwcalc(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ae_int_t n1,
     ae_int_t n2,
     ae_state *_state);
static void ftbase_internalcomplexlintranspose(/* Real    */ ae_vector* a,
     ae_int_t m,
     ae_int_t n,
     ae_int_t astart,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
static void ftbase_ffticltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
static void ftbase_fftirltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
static void ftbase_ftbasefindsmoothrec(ae_int_t n,
     ae_int_t seed,
     ae_int_t leastfactor,
     ae_int_t* best,
     ae_state *_state);









/*************************************************************************
This function is used to set error flags  during  unit  tests.  When  COND
parameter is True, FLAG variable is  set  to  True.  When  COND is  False,
FLAG is unchanged.

The purpose of this function is to have single  point  where  failures  of
unit tests can be detected.

This function returns value of COND.
*************************************************************************/
ae_bool seterrorflag(ae_bool* flag, ae_bool cond, ae_state *_state)
{
    ae_bool result;


    if( cond )
    {
        *flag = ae_true;
    }
    result = cond;
    return result;
}


/*************************************************************************
Internally calls SetErrorFlag() with condition:

    Abs(Val-RefVal)>Tol*Max(Abs(RefVal),S)
    
This function is used to test relative error in Val against  RefVal,  with
relative error being replaced by absolute when scale  of  RefVal  is  less
than S.

This function returns value of COND.
*************************************************************************/
ae_bool seterrorflagdiff(ae_bool* flag,
     double val,
     double refval,
     double tol,
     double s,
     ae_state *_state)
{
    ae_bool result;


    result = seterrorflag(flag, ae_fp_greater(ae_fabs(val-refval, _state),tol*ae_maxreal(ae_fabs(refval, _state), s, _state)), _state);
    return result;
}


/*************************************************************************
The function "touches" integer - it is used  to  avoid  compiler  messages
about unused variables (in rare cases when we do NOT want to remove  these
variables).

  -- ALGLIB --
     Copyright 17.09.2012 by Bochkanov Sergey
*************************************************************************/
void touchint(ae_int_t* a, ae_state *_state)
{


}


/*************************************************************************
The function "touches" real   -  it is used  to  avoid  compiler  messages
about unused variables (in rare cases when we do NOT want to remove  these
variables).

  -- ALGLIB --
     Copyright 17.09.2012 by Bochkanov Sergey
*************************************************************************/
void touchreal(double* a, ae_state *_state)
{


}


/*************************************************************************
The function performs zero-coalescing on real value.

NOTE: no check is performed for B<>0

  -- ALGLIB --
     Copyright 18.05.2015 by Bochkanov Sergey
*************************************************************************/
double coalesce(double a, double b, ae_state *_state)
{
    double result;


    result = a;
    if( ae_fp_eq(a,0.0) )
    {
        result = b;
    }
    return result;
}


/*************************************************************************
The function convert integer value to real value.

  -- ALGLIB --
     Copyright 17.09.2012 by Bochkanov Sergey
*************************************************************************/
double inttoreal(ae_int_t a, ae_state *_state)
{
    double result;


    result = (double)(a);
    return result;
}


/*************************************************************************
The function calculates binary logarithm.

NOTE: it costs twice as much as Ln(x)

  -- ALGLIB --
     Copyright 17.09.2012 by Bochkanov Sergey
*************************************************************************/
double logbase2(double x, ae_state *_state)
{
    double result;


    result = ae_log(x, _state)/ae_log((double)(2), _state);
    return result;
}


/*************************************************************************
This function compares two numbers for approximate equality, with tolerance
to errors as large as max(|a|,|b|)*tol.


  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool approxequalrel(double a, double b, double tol, ae_state *_state)
{
    ae_bool result;


    result = ae_fp_less_eq(ae_fabs(a-b, _state),ae_maxreal(ae_fabs(a, _state), ae_fabs(b, _state), _state)*tol);
    return result;
}


/*************************************************************************
This  function  generates  1-dimensional  general  interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1d(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    double h;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolationEqdist1D: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        x->ptr.p_double[0] = a;
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
        h = (b-a)/(n-1);
        for(i=1; i<=n-1; i++)
        {
            if( i!=n-1 )
            {
                x->ptr.p_double[i] = a+(i+0.2*(2*ae_randomreal(_state)-1))*h;
            }
            else
            {
                x->ptr.p_double[i] = b;
            }
            y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function generates  1-dimensional equidistant interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1dequidist(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    double h;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolationEqdist1D: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        x->ptr.p_double[0] = a;
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
        h = (b-a)/(n-1);
        for(i=1; i<=n-1; i++)
        {
            x->ptr.p_double[i] = a+i*h;
            y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*h;
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function generates  1-dimensional Chebyshev-1 interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1dcheb1(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolation1DCheb1: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        for(i=0; i<=n-1; i++)
        {
            x->ptr.p_double[i] = 0.5*(b+a)+0.5*(b-a)*ae_cos(ae_pi*(2*i+1)/(2*n), _state);
            if( i==0 )
            {
                y->ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            else
            {
                y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
            }
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function generates  1-dimensional Chebyshev-2 interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1dcheb2(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolation1DCheb2: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        for(i=0; i<=n-1; i++)
        {
            x->ptr.p_double[i] = 0.5*(b+a)+0.5*(b-a)*ae_cos(ae_pi*i/(n-1), _state);
            if( i==0 )
            {
                y->ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            else
            {
                y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
            }
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function checks that all values from X[] are distinct. It does more
than just usual floating point comparison:
* first, it calculates max(X) and min(X)
* second, it maps X[] from [min,max] to [1,2]
* only at this stage actual comparison is done

The meaning of such check is to ensure that all values are "distinct enough"
and will not cause interpolation subroutine to fail.

NOTE:
    X[] must be sorted by ascending (subroutine ASSERT's it)

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool aredistinct(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{
    double a;
    double b;
    ae_int_t i;
    ae_bool nonsorted;
    ae_bool result;


    ae_assert(n>=1, "APSERVAreDistinct: internal error (N<1)", _state);
    if( n==1 )
    {
        
        /*
         * everything is alright, it is up to caller to decide whether it
         * can interpolate something with just one point
         */
        result = ae_true;
        return result;
    }
    a = x->ptr.p_double[0];
    b = x->ptr.p_double[0];
    nonsorted = ae_false;
    for(i=1; i<=n-1; i++)
    {
        a = ae_minreal(a, x->ptr.p_double[i], _state);
        b = ae_maxreal(b, x->ptr.p_double[i], _state);
        nonsorted = nonsorted||ae_fp_greater_eq(x->ptr.p_double[i-1],x->ptr.p_double[i]);
    }
    ae_assert(!nonsorted, "APSERVAreDistinct: internal error (not sorted)", _state);
    for(i=1; i<=n-1; i++)
    {
        if( ae_fp_eq((x->ptr.p_double[i]-a)/(b-a)+1,(x->ptr.p_double[i-1]-a)/(b-a)+1) )
        {
            result = ae_false;
            return result;
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that two boolean values are the same (both  are  True 
or both are False).

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool aresameboolean(ae_bool v1, ae_bool v2, ae_state *_state)
{
    ae_bool result;


    result = (v1&&v2)||(!v1&&!v2);
    return result;
}


/*************************************************************************
If Length(X)<N, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void bvectorsetlengthatleast(/* Boolean */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{


    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
}


/*************************************************************************
If Length(X)<N, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void ivectorsetlengthatleast(/* Integer */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{


    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
}


/*************************************************************************
If Length(X)<N, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void rvectorsetlengthatleast(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{


    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
}


/*************************************************************************
If Cols(X)<N or Rows(X)<M, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void rmatrixsetlengthatleast(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{


    if( m>0&&n>0 )
    {
        if( x->rows<m||x->cols<n )
        {
            ae_matrix_set_length(x, m, n, _state);
        }
    }
}


/*************************************************************************
Resizes X and:
* preserves old contents of X
* fills new elements by zeros

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void rmatrixresize(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix oldx;
    ae_int_t i;
    ae_int_t j;
    ae_int_t m2;
    ae_int_t n2;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&oldx, 0, 0, DT_REAL, _state);

    m2 = x->rows;
    n2 = x->cols;
    ae_swap_matrices(x, &oldx);
    ae_matrix_set_length(x, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( i<m2&&j<n2 )
            {
                x->ptr.pp_double[i][j] = oldx.ptr.pp_double[i][j];
            }
            else
            {
                x->ptr.pp_double[i][j] = 0.0;
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Resizes X and:
* preserves old contents of X
* fills new elements by zeros

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void imatrixresize(/* Integer */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix oldx;
    ae_int_t i;
    ae_int_t j;
    ae_int_t m2;
    ae_int_t n2;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&oldx, 0, 0, DT_INT, _state);

    m2 = x->rows;
    n2 = x->cols;
    ae_swap_matrices(x, &oldx);
    ae_matrix_set_length(x, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( i<m2&&j<n2 )
            {
                x->ptr.pp_int[i][j] = oldx.ptr.pp_int[i][j];
            }
            else
            {
                x->ptr.pp_int[i][j] = 0;
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This function checks that length(X) is at least N and first N values  from
X[] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitevector(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteVector: internal error (N<0)", _state);
    if( n==0 )
    {
        result = ae_true;
        return result;
    }
    if( x->cnt<n )
    {
        result = ae_false;
        return result;
    }
    for(i=0; i<=n-1; i++)
    {
        if( !ae_isfinite(x->ptr.p_double[i], _state) )
        {
            result = ae_false;
            return result;
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that first N values from X[] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitecvector(/* Complex */ ae_vector* z,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteCVector: internal error (N<0)", _state);
    for(i=0; i<=n-1; i++)
    {
        if( !ae_isfinite(z->ptr.p_complex[i].x, _state)||!ae_isfinite(z->ptr.p_complex[i].y, _state) )
        {
            result = ae_false;
            return result;
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that size of X is at least MxN and values from
X[0..M-1,0..N-1] are finite.

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitematrix(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteMatrix: internal error (N<0)", _state);
    ae_assert(m>=0, "APSERVIsFiniteMatrix: internal error (M<0)", _state);
    if( m==0||n==0 )
    {
        result = ae_true;
        return result;
    }
    if( x->rows<m||x->cols<n )
    {
        result = ae_false;
        return result;
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( !ae_isfinite(x->ptr.pp_double[i][j], _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from X[0..M-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitecmatrix(/* Complex */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteCMatrix: internal error (N<0)", _state);
    ae_assert(m>=0, "APSERVIsFiniteCMatrix: internal error (M<0)", _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( !ae_isfinite(x->ptr.pp_complex[i][j].x, _state)||!ae_isfinite(x->ptr.pp_complex[i][j].y, _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that size of X is at least NxN and all values from
upper/lower triangle of X[0..N-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitertrmatrix(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j1;
    ae_int_t j2;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteRTRMatrix: internal error (N<0)", _state);
    if( n==0 )
    {
        result = ae_true;
        return result;
    }
    if( x->rows<n||x->cols<n )
    {
        result = ae_false;
        return result;
    }
    for(i=0; i<=n-1; i++)
    {
        if( isupper )
        {
            j1 = i;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i;
        }
        for(j=j1; j<=j2; j++)
        {
            if( !ae_isfinite(x->ptr.pp_double[i][j], _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from upper/lower triangle of
X[0..N-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitectrmatrix(/* Complex */ ae_matrix* x,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j1;
    ae_int_t j2;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteCTRMatrix: internal error (N<0)", _state);
    for(i=0; i<=n-1; i++)
    {
        if( isupper )
        {
            j1 = i;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i;
        }
        for(j=j1; j<=j2; j++)
        {
            if( !ae_isfinite(x->ptr.pp_complex[i][j].x, _state)||!ae_isfinite(x->ptr.pp_complex[i][j].y, _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from X[0..M-1,0..N-1] are  finite  or
NaN's.

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfiniteornanmatrix(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteOrNaNMatrix: internal error (N<0)", _state);
    ae_assert(m>=0, "APSERVIsFiniteOrNaNMatrix: internal error (M<0)", _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( !(ae_isfinite(x->ptr.pp_double[i][j], _state)||ae_isnan(x->ptr.pp_double[i][j], _state)) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
Safe sqrt(x^2+y^2)

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safepythag2(double x, double y, ae_state *_state)
{
    double w;
    double xabs;
    double yabs;
    double z;
    double result;


    xabs = ae_fabs(x, _state);
    yabs = ae_fabs(y, _state);
    w = ae_maxreal(xabs, yabs, _state);
    z = ae_minreal(xabs, yabs, _state);
    if( ae_fp_eq(z,(double)(0)) )
    {
        result = w;
    }
    else
    {
        result = w*ae_sqrt(1+ae_sqr(z/w, _state), _state);
    }
    return result;
}


/*************************************************************************
Safe sqrt(x^2+y^2)

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safepythag3(double x, double y, double z, ae_state *_state)
{
    double w;
    double result;


    w = ae_maxreal(ae_fabs(x, _state), ae_maxreal(ae_fabs(y, _state), ae_fabs(z, _state), _state), _state);
    if( ae_fp_eq(w,(double)(0)) )
    {
        result = (double)(0);
        return result;
    }
    x = x/w;
    y = y/w;
    z = z/w;
    result = w*ae_sqrt(ae_sqr(x, _state)+ae_sqr(y, _state)+ae_sqr(z, _state), _state);
    return result;
}


/*************************************************************************
Safe division.

This function attempts to calculate R=X/Y without overflow.

It returns:
* +1, if abs(X/Y)>=MaxRealNumber or undefined - overflow-like situation
      (no overlfow is generated, R is either NAN, PosINF, NegINF)
*  0, if MinRealNumber<abs(X/Y)<MaxRealNumber or X=0, Y<>0
      (R contains result, may be zero)
* -1, if 0<abs(X/Y)<MinRealNumber - underflow-like situation
      (R contains zero; it corresponds to underflow)

No overflow is generated in any case.

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
ae_int_t saferdiv(double x, double y, double* r, ae_state *_state)
{
    ae_int_t result;

    *r = 0;

    
    /*
     * Two special cases:
     * * Y=0
     * * X=0 and Y<>0
     */
    if( ae_fp_eq(y,(double)(0)) )
    {
        result = 1;
        if( ae_fp_eq(x,(double)(0)) )
        {
            *r = _state->v_nan;
        }
        if( ae_fp_greater(x,(double)(0)) )
        {
            *r = _state->v_posinf;
        }
        if( ae_fp_less(x,(double)(0)) )
        {
            *r = _state->v_neginf;
        }
        return result;
    }
    if( ae_fp_eq(x,(double)(0)) )
    {
        *r = (double)(0);
        result = 0;
        return result;
    }
    
    /*
     * make Y>0
     */
    if( ae_fp_less(y,(double)(0)) )
    {
        x = -x;
        y = -y;
    }
    
    /*
     *
     */
    if( ae_fp_greater_eq(y,(double)(1)) )
    {
        *r = x/y;
        if( ae_fp_less_eq(ae_fabs(*r, _state),ae_minrealnumber) )
        {
            result = -1;
            *r = (double)(0);
        }
        else
        {
            result = 0;
        }
    }
    else
    {
        if( ae_fp_greater_eq(ae_fabs(x, _state),ae_maxrealnumber*y) )
        {
            if( ae_fp_greater(x,(double)(0)) )
            {
                *r = _state->v_posinf;
            }
            else
            {
                *r = _state->v_neginf;
            }
            result = 1;
        }
        else
        {
            *r = x/y;
            result = 0;
        }
    }
    return result;
}


/*************************************************************************
This function calculates "safe" min(X/Y,V) for positive finite X, Y, V.
No overflow is generated in any case.

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safeminposrv(double x, double y, double v, ae_state *_state)
{
    double r;
    double result;


    if( ae_fp_greater_eq(y,(double)(1)) )
    {
        
        /*
         * Y>=1, we can safely divide by Y
         */
        r = x/y;
        result = v;
        if( ae_fp_greater(v,r) )
        {
            result = r;
        }
        else
        {
            result = v;
        }
    }
    else
    {
        
        /*
         * Y<1, we can safely multiply by Y
         */
        if( ae_fp_less(x,v*y) )
        {
            result = x/y;
        }
        else
        {
            result = v;
        }
    }
    return result;
}


/*************************************************************************
This function makes periodic mapping of X to [A,B].

It accepts X, A, B (A>B). It returns T which lies in  [A,B] and integer K,
such that X = T + K*(B-A).

NOTES:
* K is represented as real value, although actually it is integer
* T is guaranteed to be in [A,B]
* T replaces X

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
void apperiodicmap(double* x,
     double a,
     double b,
     double* k,
     ae_state *_state)
{

    *k = 0;

    ae_assert(ae_fp_less(a,b), "APPeriodicMap: internal error!", _state);
    *k = (double)(ae_ifloor((*x-a)/(b-a), _state));
    *x = *x-*k*(b-a);
    while(ae_fp_less(*x,a))
    {
        *x = *x+(b-a);
        *k = *k-1;
    }
    while(ae_fp_greater(*x,b))
    {
        *x = *x-(b-a);
        *k = *k+1;
    }
    *x = ae_maxreal(*x, a, _state);
    *x = ae_minreal(*x, b, _state);
}


/*************************************************************************
Returns random normal number using low-quality system-provided generator

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
double randomnormal(ae_state *_state)
{
    double u;
    double v;
    double s;
    double result;


    for(;;)
    {
        u = 2*ae_randomreal(_state)-1;
        v = 2*ae_randomreal(_state)-1;
        s = ae_sqr(u, _state)+ae_sqr(v, _state);
        if( ae_fp_greater(s,(double)(0))&&ae_fp_less(s,(double)(1)) )
        {
            
            /*
             * two Sqrt's instead of one to
             * avoid overflow when S is too small
             */
            s = ae_sqrt(-2*ae_log(s, _state), _state)/ae_sqrt(s, _state);
            result = u*s;
            break;
        }
    }
    return result;
}


/*************************************************************************
Generates random unit vector using low-quality system-provided generator.
Reallocates array if its size is too short.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void randomunit(ae_int_t n, /* Real    */ ae_vector* x, ae_state *_state)
{
    ae_int_t i;
    double v;
    double vv;


    ae_assert(n>0, "RandomUnit: N<=0", _state);
    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
    do
    {
        v = 0.0;
        for(i=0; i<=n-1; i++)
        {
            vv = randomnormal(_state);
            x->ptr.p_double[i] = vv;
            v = v+vv*vv;
        }
    }
    while(ae_fp_less_eq(v,(double)(0)));
    v = 1/ae_sqrt(v, _state);
    for(i=0; i<=n-1; i++)
    {
        x->ptr.p_double[i] = x->ptr.p_double[i]*v;
    }
}


/*************************************************************************
This function is used to swap two integer values
*************************************************************************/
void swapi(ae_int_t* v0, ae_int_t* v1, ae_state *_state)
{
    ae_int_t v;


    v = *v0;
    *v0 = *v1;
    *v1 = v;
}


/*************************************************************************
This function is used to swap two real values
*************************************************************************/
void swapr(double* v0, double* v1, ae_state *_state)
{
    double v;


    v = *v0;
    *v0 = *v1;
    *v1 = v;
}


/*************************************************************************
This function is used to return maximum of three real values
*************************************************************************/
double maxreal3(double v0, double v1, double v2, ae_state *_state)
{
    double result;


    result = v0;
    if( ae_fp_less(result,v1) )
    {
        result = v1;
    }
    if( ae_fp_less(result,v2) )
    {
        result = v2;
    }
    return result;
}


/*************************************************************************
This function is used to increment value of integer variable
*************************************************************************/
void inc(ae_int_t* v, ae_state *_state)
{


    *v = *v+1;
}


/*************************************************************************
This function is used to decrement value of integer variable
*************************************************************************/
void dec(ae_int_t* v, ae_state *_state)
{


    *v = *v-1;
}


/*************************************************************************
This function performs two operations:
1. decrements value of integer variable, if it is positive
2. explicitly sets variable to zero if it is non-positive
It is used by some algorithms to decrease value of internal counters.
*************************************************************************/
void countdown(ae_int_t* v, ae_state *_state)
{


    if( *v>0 )
    {
        *v = *v-1;
    }
    else
    {
        *v = 0;
    }
}


/*************************************************************************
'bounds' value: maps X to [B1,B2]

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
double boundval(double x, double b1, double b2, ae_state *_state)
{
    double result;


    if( ae_fp_less_eq(x,b1) )
    {
        result = b1;
        return result;
    }
    if( ae_fp_greater_eq(x,b2) )
    {
        result = b2;
        return result;
    }
    result = x;
    return result;
}


/*************************************************************************
Allocation of serializer: complex value
*************************************************************************/
void alloccomplex(ae_serializer* s, ae_complex v, ae_state *_state)
{


    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
}


/*************************************************************************
Serialization: complex value
*************************************************************************/
void serializecomplex(ae_serializer* s, ae_complex v, ae_state *_state)
{


    ae_serializer_serialize_double(s, v.x, _state);
    ae_serializer_serialize_double(s, v.y, _state);
}


/*************************************************************************
Unserialization: complex value
*************************************************************************/
ae_complex unserializecomplex(ae_serializer* s, ae_state *_state)
{
    ae_complex result;


    ae_serializer_unserialize_double(s, &result.x, _state);
    ae_serializer_unserialize_double(s, &result.y, _state);
    return result;
}


/*************************************************************************
Allocation of serializer: real array
*************************************************************************/
void allocrealarray(ae_serializer* s,
     /* Real    */ ae_vector* v,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;


    if( n<0 )
    {
        n = v->cnt;
    }
    ae_serializer_alloc_entry(s);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_alloc_entry(s);
    }
}


/*************************************************************************
Serialization: complex value
*************************************************************************/
void serializerealarray(ae_serializer* s,
     /* Real    */ ae_vector* v,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;


    if( n<0 )
    {
        n = v->cnt;
    }
    ae_serializer_serialize_int(s, n, _state);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_serialize_double(s, v->ptr.p_double[i], _state);
    }
}


/*************************************************************************
Unserialization: complex value
*************************************************************************/
void unserializerealarray(ae_serializer* s,
     /* Real    */ ae_vector* v,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    double t;

    ae_vector_clear(v);

    ae_serializer_unserialize_int(s, &n, _state);
    if( n==0 )
    {
        return;
    }
    ae_vector_set_length(v, n, _state);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_unserialize_double(s, &t, _state);
        v->ptr.p_double[i] = t;
    }
}


/*************************************************************************
Allocation of serializer: Integer array
*************************************************************************/
void allocintegerarray(ae_serializer* s,
     /* Integer */ ae_vector* v,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;


    if( n<0 )
    {
        n = v->cnt;
    }
    ae_serializer_alloc_entry(s);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_alloc_entry(s);
    }
}


/*************************************************************************
Serialization: Integer array
*************************************************************************/
void serializeintegerarray(ae_serializer* s,
     /* Integer */ ae_vector* v,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;


    if( n<0 )
    {
        n = v->cnt;
    }
    ae_serializer_serialize_int(s, n, _state);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_serialize_int(s, v->ptr.p_int[i], _state);
    }
}


/*************************************************************************
Unserialization: complex value
*************************************************************************/
void unserializeintegerarray(ae_serializer* s,
     /* Integer */ ae_vector* v,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    ae_int_t t;

    ae_vector_clear(v);

    ae_serializer_unserialize_int(s, &n, _state);
    if( n==0 )
    {
        return;
    }
    ae_vector_set_length(v, n, _state);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_unserialize_int(s, &t, _state);
        v->ptr.p_int[i] = t;
    }
}


/*************************************************************************
Allocation of serializer: real matrix
*************************************************************************/
void allocrealmatrix(ae_serializer* s,
     /* Real    */ ae_matrix* v,
     ae_int_t n0,
     ae_int_t n1,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    if( n0<0 )
    {
        n0 = v->rows;
    }
    if( n1<0 )
    {
        n1 = v->cols;
    }
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    for(i=0; i<=n0-1; i++)
    {
        for(j=0; j<=n1-1; j++)
        {
            ae_serializer_alloc_entry(s);
        }
    }
}


/*************************************************************************
Serialization: complex value
*************************************************************************/
void serializerealmatrix(ae_serializer* s,
     /* Real    */ ae_matrix* v,
     ae_int_t n0,
     ae_int_t n1,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    if( n0<0 )
    {
        n0 = v->rows;
    }
    if( n1<0 )
    {
        n1 = v->cols;
    }
    ae_serializer_serialize_int(s, n0, _state);
    ae_serializer_serialize_int(s, n1, _state);
    for(i=0; i<=n0-1; i++)
    {
        for(j=0; j<=n1-1; j++)
        {
            ae_serializer_serialize_double(s, v->ptr.pp_double[i][j], _state);
        }
    }
}


/*************************************************************************
Unserialization: complex value
*************************************************************************/
void unserializerealmatrix(ae_serializer* s,
     /* Real    */ ae_matrix* v,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t n0;
    ae_int_t n1;
    double t;

    ae_matrix_clear(v);

    ae_serializer_unserialize_int(s, &n0, _state);
    ae_serializer_unserialize_int(s, &n1, _state);
    if( n0==0||n1==0 )
    {
        return;
    }
    ae_matrix_set_length(v, n0, n1, _state);
    for(i=0; i<=n0-1; i++)
    {
        for(j=0; j<=n1-1; j++)
        {
            ae_serializer_unserialize_double(s, &t, _state);
            v->ptr.pp_double[i][j] = t;
        }
    }
}


/*************************************************************************
Copy integer array
*************************************************************************/
void copyintegerarray(/* Integer */ ae_vector* src,
     /* Integer */ ae_vector* dst,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(dst);

    if( src->cnt>0 )
    {
        ae_vector_set_length(dst, src->cnt, _state);
        for(i=0; i<=src->cnt-1; i++)
        {
            dst->ptr.p_int[i] = src->ptr.p_int[i];
        }
    }
}


/*************************************************************************
Copy real array
*************************************************************************/
void copyrealarray(/* Real    */ ae_vector* src,
     /* Real    */ ae_vector* dst,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(dst);

    if( src->cnt>0 )
    {
        ae_vector_set_length(dst, src->cnt, _state);
        for(i=0; i<=src->cnt-1; i++)
        {
            dst->ptr.p_double[i] = src->ptr.p_double[i];
        }
    }
}


/*************************************************************************
Copy real matrix
*************************************************************************/
void copyrealmatrix(/* Real    */ ae_matrix* src,
     /* Real    */ ae_matrix* dst,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;

    ae_matrix_clear(dst);

    if( src->rows>0&&src->cols>0 )
    {
        ae_matrix_set_length(dst, src->rows, src->cols, _state);
        for(i=0; i<=src->rows-1; i++)
        {
            for(j=0; j<=src->cols-1; j++)
            {
                dst->ptr.pp_double[i][j] = src->ptr.pp_double[i][j];
            }
        }
    }
}


/*************************************************************************
This function searches integer array. Elements in this array are actually
records, each NRec elements wide. Each record has unique header - NHeader
integer values, which identify it. Records are lexicographically sorted by
header.

Records are identified by their index, not offset (offset = NRec*index).

This function searches A (records with indices [I0,I1)) for a record with
header B. It returns index of this record (not offset!), or -1 on failure.

  -- ALGLIB --
     Copyright 28.03.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t recsearch(/* Integer */ ae_vector* a,
     ae_int_t nrec,
     ae_int_t nheader,
     ae_int_t i0,
     ae_int_t i1,
     /* Integer */ ae_vector* b,
     ae_state *_state)
{
    ae_int_t mididx;
    ae_int_t cflag;
    ae_int_t k;
    ae_int_t offs;
    ae_int_t result;


    result = -1;
    for(;;)
    {
        if( i0>=i1 )
        {
            break;
        }
        mididx = (i0+i1)/2;
        offs = nrec*mididx;
        cflag = 0;
        for(k=0; k<=nheader-1; k++)
        {
            if( a->ptr.p_int[offs+k]<b->ptr.p_int[k] )
            {
                cflag = -1;
                break;
            }
            if( a->ptr.p_int[offs+k]>b->ptr.p_int[k] )
            {
                cflag = 1;
                break;
            }
        }
        if( cflag==0 )
        {
            result = mididx;
            return result;
        }
        if( cflag<0 )
        {
            i0 = mididx+1;
        }
        else
        {
            i1 = mididx;
        }
    }
    return result;
}


/*************************************************************************
This function is used in parallel functions for recurrent division of large
task into two smaller tasks.

It has following properties:
* it works only for TaskSize>=2 (assertion is thrown otherwise)
* for TaskSize=2, it returns Task0=1, Task1=1
* in case TaskSize is odd,  Task0=TaskSize-1, Task1=1
* in case TaskSize is even, Task0 and Task1 are approximately TaskSize/2
  and both Task0 and Task1 are even, Task0>=Task1

  -- ALGLIB --
     Copyright 07.04.2013 by Bochkanov Sergey
*************************************************************************/
void splitlengtheven(ae_int_t tasksize,
     ae_int_t* task0,
     ae_int_t* task1,
     ae_state *_state)
{

    *task0 = 0;
    *task1 = 0;

    ae_assert(tasksize>=2, "SplitLengthEven: TaskSize<2", _state);
    if( tasksize==2 )
    {
        *task0 = 1;
        *task1 = 1;
        return;
    }
    if( tasksize%2==0 )
    {
        
        /*
         * Even division
         */
        *task0 = tasksize/2;
        *task1 = tasksize/2;
        if( *task0%2!=0 )
        {
            *task0 = *task0+1;
            *task1 = *task1-1;
        }
    }
    else
    {
        
        /*
         * Odd task size, split trailing odd part from it.
         */
        *task0 = tasksize-1;
        *task1 = 1;
    }
    ae_assert(*task0>=1, "SplitLengthEven: internal error", _state);
    ae_assert(*task1>=1, "SplitLengthEven: internal error", _state);
}


/*************************************************************************
This function is used in parallel functions for recurrent division of large
task into two smaller tasks.

It has following properties:
* it works only for TaskSize>=2 and ChunkSize>=2
  (assertion is thrown otherwise)
* Task0+Task1=TaskSize, Task0>0, Task1>0
* Task0 and Task1 are close to each other
* in case TaskSize>ChunkSize, Task0 is always divisible by ChunkSize

  -- ALGLIB --
     Copyright 07.04.2013 by Bochkanov Sergey
*************************************************************************/
void splitlength(ae_int_t tasksize,
     ae_int_t chunksize,
     ae_int_t* task0,
     ae_int_t* task1,
     ae_state *_state)
{

    *task0 = 0;
    *task1 = 0;

    ae_assert(chunksize>=2, "SplitLength: ChunkSize<2", _state);
    ae_assert(tasksize>=2, "SplitLength: TaskSize<2", _state);
    *task0 = tasksize/2;
    if( *task0>chunksize&&*task0%chunksize!=0 )
    {
        *task0 = *task0-*task0%chunksize;
    }
    *task1 = tasksize-(*task0);
    ae_assert(*task0>=1, "SplitLength: internal error", _state);
    ae_assert(*task1>=1, "SplitLength: internal error", _state);
}


/*************************************************************************
This function is used to calculate number of chunks (including partial,
non-complete chunks) in some set. It expects that ChunkSize>=1, TaskSize>=0.
Assertion is thrown otherwise.

Function result is equivalent to Ceil(TaskSize/ChunkSize), but with guarantees
that rounding errors won't ruin results.

  -- ALGLIB --
     Copyright 21.01.2015 by Bochkanov Sergey
*************************************************************************/
ae_int_t chunkscount(ae_int_t tasksize,
     ae_int_t chunksize,
     ae_state *_state)
{
    ae_int_t result;


    ae_assert(tasksize>=0, "ChunksCount: TaskSize<0", _state);
    ae_assert(chunksize>=1, "ChunksCount: ChunkSize<1", _state);
    result = tasksize/chunksize;
    if( tasksize%chunksize!=0 )
    {
        result = result+1;
    }
    return result;
}


void _apbuffers_init(void* _p, ae_state *_state)
{
    apbuffers *p = (apbuffers*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_init(&p->ba0, 0, DT_BOOL, _state);
    ae_vector_init(&p->ia0, 0, DT_INT, _state);
    ae_vector_init(&p->ia1, 0, DT_INT, _state);
    ae_vector_init(&p->ia2, 0, DT_INT, _state);
    ae_vector_init(&p->ia3, 0, DT_INT, _state);
    ae_vector_init(&p->ra0, 0, DT_REAL, _state);
    ae_vector_init(&p->ra1, 0, DT_REAL, _state);
    ae_vector_init(&p->ra2, 0, DT_REAL, _state);
    ae_vector_init(&p->ra3, 0, DT_REAL, _state);
    ae_matrix_init(&p->rm0, 0, 0, DT_REAL, _state);
    ae_matrix_init(&p->rm1, 0, 0, DT_REAL, _state);
}


void _apbuffers_init_copy(void* _dst, void* _src, ae_state *_state)
{
    apbuffers *dst = (apbuffers*)_dst;
    apbuffers *src = (apbuffers*)_src;
    ae_vector_init_copy(&dst->ba0, &src->ba0, _state);
    ae_vector_init_copy(&dst->ia0, &src->ia0, _state);
    ae_vector_init_copy(&dst->ia1, &src->ia1, _state);
    ae_vector_init_copy(&dst->ia2, &src->ia2, _state);
    ae_vector_init_copy(&dst->ia3, &src->ia3, _state);
    ae_vector_init_copy(&dst->ra0, &src->ra0, _state);
    ae_vector_init_copy(&dst->ra1, &src->ra1, _state);
    ae_vector_init_copy(&dst->ra2, &src->ra2, _state);
    ae_vector_init_copy(&dst->ra3, &src->ra3, _state);
    ae_matrix_init_copy(&dst->rm0, &src->rm0, _state);
    ae_matrix_init_copy(&dst->rm1, &src->rm1, _state);
}


void _apbuffers_clear(void* _p)
{
    apbuffers *p = (apbuffers*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_clear(&p->ba0);
    ae_vector_clear(&p->ia0);
    ae_vector_clear(&p->ia1);
    ae_vector_clear(&p->ia2);
    ae_vector_clear(&p->ia3);
    ae_vector_clear(&p->ra0);
    ae_vector_clear(&p->ra1);
    ae_vector_clear(&p->ra2);
    ae_vector_clear(&p->ra3);
    ae_matrix_clear(&p->rm0);
    ae_matrix_clear(&p->rm1);
}


void _apbuffers_destroy(void* _p)
{
    apbuffers *p = (apbuffers*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_destroy(&p->ba0);
    ae_vector_destroy(&p->ia0);
    ae_vector_destroy(&p->ia1);
    ae_vector_destroy(&p->ia2);
    ae_vector_destroy(&p->ia3);
    ae_vector_destroy(&p->ra0);
    ae_vector_destroy(&p->ra1);
    ae_vector_destroy(&p->ra2);
    ae_vector_destroy(&p->ra3);
    ae_matrix_destroy(&p->rm0);
    ae_matrix_destroy(&p->rm1);
}


void _sboolean_init(void* _p, ae_state *_state)
{
    sboolean *p = (sboolean*)_p;
    ae_touch_ptr((void*)p);
}


void _sboolean_init_copy(void* _dst, void* _src, ae_state *_state)
{
    sboolean *dst = (sboolean*)_dst;
    sboolean *src = (sboolean*)_src;
    dst->val = src->val;
}


void _sboolean_clear(void* _p)
{
    sboolean *p = (sboolean*)_p;
    ae_touch_ptr((void*)p);
}


void _sboolean_destroy(void* _p)
{
    sboolean *p = (sboolean*)_p;
    ae_touch_ptr((void*)p);
}


void _sbooleanarray_init(void* _p, ae_state *_state)
{
    sbooleanarray *p = (sbooleanarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_init(&p->val, 0, DT_BOOL, _state);
}


void _sbooleanarray_init_copy(void* _dst, void* _src, ae_state *_state)
{
    sbooleanarray *dst = (sbooleanarray*)_dst;
    sbooleanarray *src = (sbooleanarray*)_src;
    ae_vector_init_copy(&dst->val, &src->val, _state);
}


void _sbooleanarray_clear(void* _p)
{
    sbooleanarray *p = (sbooleanarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_clear(&p->val);
}


void _sbooleanarray_destroy(void* _p)
{
    sbooleanarray *p = (sbooleanarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_destroy(&p->val);
}


void _sinteger_init(void* _p, ae_state *_state)
{
    sinteger *p = (sinteger*)_p;
    ae_touch_ptr((void*)p);
}


void _sinteger_init_copy(void* _dst, void* _src, ae_state *_state)
{
    sinteger *dst = (sinteger*)_dst;
    sinteger *src = (sinteger*)_src;
    dst->val = src->val;
}


void _sinteger_clear(void* _p)
{
    sinteger *p = (sinteger*)_p;
    ae_touch_ptr((void*)p);
}


void _sinteger_destroy(void* _p)
{
    sinteger *p = (sinteger*)_p;
    ae_touch_ptr((void*)p);
}


void _sintegerarray_init(void* _p, ae_state *_state)
{
    sintegerarray *p = (sintegerarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_init(&p->val, 0, DT_INT, _state);
}


void _sintegerarray_init_copy(void* _dst, void* _src, ae_state *_state)
{
    sintegerarray *dst = (sintegerarray*)_dst;
    sintegerarray *src = (sintegerarray*)_src;
    ae_vector_init_copy(&dst->val, &src->val, _state);
}


void _sintegerarray_clear(void* _p)
{
    sintegerarray *p = (sintegerarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_clear(&p->val);
}


void _sintegerarray_destroy(void* _p)
{
    sintegerarray *p = (sintegerarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_destroy(&p->val);
}


void _sreal_init(void* _p, ae_state *_state)
{
    sreal *p = (sreal*)_p;
    ae_touch_ptr((void*)p);
}


void _sreal_init_copy(void* _dst, void* _src, ae_state *_state)
{
    sreal *dst = (sreal*)_dst;
    sreal *src = (sreal*)_src;
    dst->val = src->val;
}


void _sreal_clear(void* _p)
{
    sreal *p = (sreal*)_p;
    ae_touch_ptr((void*)p);
}


void _sreal_destroy(void* _p)
{
    sreal *p = (sreal*)_p;
    ae_touch_ptr((void*)p);
}


void _srealarray_init(void* _p, ae_state *_state)
{
    srealarray *p = (srealarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_init(&p->val, 0, DT_REAL, _state);
}


void _srealarray_init_copy(void* _dst, void* _src, ae_state *_state)
{
    srealarray *dst = (srealarray*)_dst;
    srealarray *src = (srealarray*)_src;
    ae_vector_init_copy(&dst->val, &src->val, _state);
}


void _srealarray_clear(void* _p)
{
    srealarray *p = (srealarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_clear(&p->val);
}


void _srealarray_destroy(void* _p)
{
    srealarray *p = (srealarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_destroy(&p->val);
}


void _scomplex_init(void* _p, ae_state *_state)
{
    scomplex *p = (scomplex*)_p;
    ae_touch_ptr((void*)p);
}


void _scomplex_init_copy(void* _dst, void* _src, ae_state *_state)
{
    scomplex *dst = (scomplex*)_dst;
    scomplex *src = (scomplex*)_src;
    dst->val = src->val;
}


void _scomplex_clear(void* _p)
{
    scomplex *p = (scomplex*)_p;
    ae_touch_ptr((void*)p);
}


void _scomplex_destroy(void* _p)
{
    scomplex *p = (scomplex*)_p;
    ae_touch_ptr((void*)p);
}


void _scomplexarray_init(void* _p, ae_state *_state)
{
    scomplexarray *p = (scomplexarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_init(&p->val, 0, DT_COMPLEX, _state);
}


void _scomplexarray_init_copy(void* _dst, void* _src, ae_state *_state)
{
    scomplexarray *dst = (scomplexarray*)_dst;
    scomplexarray *src = (scomplexarray*)_src;
    ae_vector_init_copy(&dst->val, &src->val, _state);
}


void _scomplexarray_clear(void* _p)
{
    scomplexarray *p = (scomplexarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_clear(&p->val);
}


void _scomplexarray_destroy(void* _p)
{
    scomplexarray *p = (scomplexarray*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_destroy(&p->val);
}




ae_int_t getrdfserializationcode(ae_state *_state)
{
    ae_int_t result;


    result = 1;
    return result;
}


ae_int_t getkdtreeserializationcode(ae_state *_state)
{
    ae_int_t result;


    result = 2;
    return result;
}


ae_int_t getmlpserializationcode(ae_state *_state)
{
    ae_int_t result;


    result = 3;
    return result;
}


ae_int_t getmlpeserializationcode(ae_state *_state)
{
    ae_int_t result;


    result = 4;
    return result;
}


ae_int_t getrbfserializationcode(ae_state *_state)
{
    ae_int_t result;


    result = 5;
    return result;
}




/*************************************************************************
This function sorts array of real keys by ascending.

Its results are:
* sorted array A
* permutation tables P1, P2

Algorithm outputs permutation tables using two formats:
* as usual permutation of [0..N-1]. If P1[i]=j, then sorted A[i]  contains
  value which was moved there from J-th position.
* as a sequence of pairwise permutations. Sorted A[] may  be  obtained  by
  swaping A[i] and A[P2[i]] for all i from 0 to N-1.
  
INPUT PARAMETERS:
    A       -   unsorted array
    N       -   array size

OUPUT PARAMETERS:
    A       -   sorted array
    P1, P2  -   permutation tables, array[N]
    
NOTES:
    this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.

  -- ALGLIB --
     Copyright 14.05.2008 by Bochkanov Sergey
*************************************************************************/
void tagsort(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     ae_state *_state)
{
    ae_frame _frame_block;
    apbuffers buf;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(p1);
    ae_vector_clear(p2);
    _apbuffers_init(&buf, _state);

    tagsortbuf(a, n, p1, p2, &buf, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Buffered variant of TagSort, which accepts preallocated output arrays as
well as special structure for buffered allocations. If arrays are too
short, they are reallocated. If they are large enough, no memory
allocation is done.

It is intended to be used in the performance-critical parts of code, where
additional allocations can lead to severe performance degradation

  -- ALGLIB --
     Copyright 14.05.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortbuf(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     apbuffers* buf,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t lv;
    ae_int_t lp;
    ae_int_t rv;
    ae_int_t rp;


    
    /*
     * Special cases
     */
    if( n<=0 )
    {
        return;
    }
    if( n==1 )
    {
        ivectorsetlengthatleast(p1, 1, _state);
        ivectorsetlengthatleast(p2, 1, _state);
        p1->ptr.p_int[0] = 0;
        p2->ptr.p_int[0] = 0;
        return;
    }
    
    /*
     * General case, N>1: prepare permutations table P1
     */
    ivectorsetlengthatleast(p1, n, _state);
    for(i=0; i<=n-1; i++)
    {
        p1->ptr.p_int[i] = i;
    }
    
    /*
     * General case, N>1: sort, update P1
     */
    rvectorsetlengthatleast(&buf->ra0, n, _state);
    ivectorsetlengthatleast(&buf->ia0, n, _state);
    tagsortfasti(a, p1, &buf->ra0, &buf->ia0, n, _state);
    
    /*
     * General case, N>1: fill permutations table P2
     *
     * To fill P2 we maintain two arrays:
     * * PV (Buf.IA0), Position(Value). PV[i] contains position of I-th key at the moment
     * * VP (Buf.IA1), Value(Position). VP[i] contains key which has position I at the moment
     *
     * At each step we making permutation of two items:
     *   Left, which is given by position/value pair LP/LV
     *   and Right, which is given by RP/RV
     * and updating PV[] and VP[] correspondingly.
     */
    ivectorsetlengthatleast(&buf->ia0, n, _state);
    ivectorsetlengthatleast(&buf->ia1, n, _state);
    ivectorsetlengthatleast(p2, n, _state);
    for(i=0; i<=n-1; i++)
    {
        buf->ia0.ptr.p_int[i] = i;
        buf->ia1.ptr.p_int[i] = i;
    }
    for(i=0; i<=n-1; i++)
    {
        
        /*
         * calculate LP, LV, RP, RV
         */
        lp = i;
        lv = buf->ia1.ptr.p_int[lp];
        rv = p1->ptr.p_int[i];
        rp = buf->ia0.ptr.p_int[rv];
        
        /*
         * Fill P2
         */
        p2->ptr.p_int[i] = rp;
        
        /*
         * update PV and VP
         */
        buf->ia1.ptr.p_int[lp] = rv;
        buf->ia1.ptr.p_int[rp] = lv;
        buf->ia0.ptr.p_int[lv] = rp;
        buf->ia0.ptr.p_int[rv] = lp;
    }
}


/*************************************************************************
Same as TagSort, but optimized for real keys and integer labels.

A is sorted, and same permutations are applied to B.

NOTES:
1.  this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.
2.  this function uses two buffers, BufA and BufB, each is N elements large.
    They may be preallocated (which will save some time) or not, in which
    case function will automatically allocate memory.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortfasti(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Integer */ ae_vector* bufb,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool isascending;
    ae_bool isdescending;
    double tmpr;
    ae_int_t tmpi;


    
    /*
     * Special case
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * Test for already sorted set
     */
    isascending = ae_true;
    isdescending = ae_true;
    for(i=1; i<=n-1; i++)
    {
        isascending = isascending&&a->ptr.p_double[i]>=a->ptr.p_double[i-1];
        isdescending = isdescending&&a->ptr.p_double[i]<=a->ptr.p_double[i-1];
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i=0; i<=n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmpr = a->ptr.p_double[i];
            a->ptr.p_double[i] = a->ptr.p_double[j];
            a->ptr.p_double[j] = tmpr;
            tmpi = b->ptr.p_int[i];
            b->ptr.p_int[i] = b->ptr.p_int[j];
            b->ptr.p_int[j] = tmpi;
        }
        return;
    }
    
    /*
     * General case
     */
    if( bufa->cnt<n )
    {
        ae_vector_set_length(bufa, n, _state);
    }
    if( bufb->cnt<n )
    {
        ae_vector_set_length(bufb, n, _state);
    }
    tsort_tagsortfastirec(a, b, bufa, bufb, 0, n-1, _state);
}


/*************************************************************************
Same as TagSort, but optimized for real keys and real labels.

A is sorted, and same permutations are applied to B.

NOTES:
1.  this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.
2.  this function uses two buffers, BufA and BufB, each is N elements large.
    They may be preallocated (which will save some time) or not, in which
    case function will automatically allocate memory.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortfastr(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool isascending;
    ae_bool isdescending;
    double tmpr;


    
    /*
     * Special case
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * Test for already sorted set
     */
    isascending = ae_true;
    isdescending = ae_true;
    for(i=1; i<=n-1; i++)
    {
        isascending = isascending&&a->ptr.p_double[i]>=a->ptr.p_double[i-1];
        isdescending = isdescending&&a->ptr.p_double[i]<=a->ptr.p_double[i-1];
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i=0; i<=n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmpr = a->ptr.p_double[i];
            a->ptr.p_double[i] = a->ptr.p_double[j];
            a->ptr.p_double[j] = tmpr;
            tmpr = b->ptr.p_double[i];
            b->ptr.p_double[i] = b->ptr.p_double[j];
            b->ptr.p_double[j] = tmpr;
        }
        return;
    }
    
    /*
     * General case
     */
    if( bufa->cnt<n )
    {
        ae_vector_set_length(bufa, n, _state);
    }
    if( bufb->cnt<n )
    {
        ae_vector_set_length(bufb, n, _state);
    }
    tsort_tagsortfastrrec(a, b, bufa, bufb, 0, n-1, _state);
}


/*************************************************************************
Same as TagSort, but optimized for real keys without labels.

A is sorted, and that's all.

NOTES:
1.  this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.
2.  this function uses buffer, BufA, which is N elements large. It may be
    preallocated (which will save some time) or not, in which case
    function will automatically allocate memory.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortfast(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* bufa,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool isascending;
    ae_bool isdescending;
    double tmpr;


    
    /*
     * Special case
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * Test for already sorted set
     */
    isascending = ae_true;
    isdescending = ae_true;
    for(i=1; i<=n-1; i++)
    {
        isascending = isascending&&a->ptr.p_double[i]>=a->ptr.p_double[i-1];
        isdescending = isdescending&&a->ptr.p_double[i]<=a->ptr.p_double[i-1];
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i=0; i<=n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmpr = a->ptr.p_double[i];
            a->ptr.p_double[i] = a->ptr.p_double[j];
            a->ptr.p_double[j] = tmpr;
        }
        return;
    }
    
    /*
     * General case
     */
    if( bufa->cnt<n )
    {
        ae_vector_set_length(bufa, n, _state);
    }
    tsort_tagsortfastrec(a, bufa, 0, n-1, _state);
}


/*************************************************************************
Sorting function optimized for integer keys and real labels, can be used
to sort middle of the array

A is sorted, and same permutations are applied to B.

NOTES:
    this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortmiddleir(/* Integer */ ae_vector* a,
     /* Real    */ ae_vector* b,
     ae_int_t offset,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;
    ae_int_t t;
    ae_int_t tmp;
    double tmpr;


    
    /*
     * Special cases
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * General case, N>1: sort, update B
     */
    i = 2;
    do
    {
        t = i;
        while(t!=1)
        {
            k = t/2;
            if( a->ptr.p_int[offset+k-1]>=a->ptr.p_int[offset+t-1] )
            {
                t = 1;
            }
            else
            {
                tmp = a->ptr.p_int[offset+k-1];
                a->ptr.p_int[offset+k-1] = a->ptr.p_int[offset+t-1];
                a->ptr.p_int[offset+t-1] = tmp;
                tmpr = b->ptr.p_double[offset+k-1];
                b->ptr.p_double[offset+k-1] = b->ptr.p_double[offset+t-1];
                b->ptr.p_double[offset+t-1] = tmpr;
                t = k;
            }
        }
        i = i+1;
    }
    while(i<=n);
    i = n-1;
    do
    {
        tmp = a->ptr.p_int[offset+i];
        a->ptr.p_int[offset+i] = a->ptr.p_int[offset+0];
        a->ptr.p_int[offset+0] = tmp;
        tmpr = b->ptr.p_double[offset+i];
        b->ptr.p_double[offset+i] = b->ptr.p_double[offset+0];
        b->ptr.p_double[offset+0] = tmpr;
        t = 1;
        while(t!=0)
        {
            k = 2*t;
            if( k>i )
            {
                t = 0;
            }
            else
            {
                if( k<i )
                {
                    if( a->ptr.p_int[offset+k]>a->ptr.p_int[offset+k-1] )
                    {
                        k = k+1;
                    }
                }
                if( a->ptr.p_int[offset+t-1]>=a->ptr.p_int[offset+k-1] )
                {
                    t = 0;
                }
                else
                {
                    tmp = a->ptr.p_int[offset+k-1];
                    a->ptr.p_int[offset+k-1] = a->ptr.p_int[offset+t-1];
                    a->ptr.p_int[offset+t-1] = tmp;
                    tmpr = b->ptr.p_double[offset+k-1];
                    b->ptr.p_double[offset+k-1] = b->ptr.p_double[offset+t-1];
                    b->ptr.p_double[offset+t-1] = tmpr;
                    t = k;
                }
            }
        }
        i = i-1;
    }
    while(i>=1);
}


/*************************************************************************
Heap operations: adds element to the heap

PARAMETERS:
    A       -   heap itself, must be at least array[0..N]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap (without new element).
                updated on output
    VA      -   value of the element being added
    VB      -   value of the tag

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheappushi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t* n,
     double va,
     ae_int_t vb,
     ae_state *_state)
{
    ae_int_t j;
    ae_int_t k;
    double v;


    if( *n<0 )
    {
        return;
    }
    
    /*
     * N=0 is a special case
     */
    if( *n==0 )
    {
        a->ptr.p_double[0] = va;
        b->ptr.p_int[0] = vb;
        *n = *n+1;
        return;
    }
    
    /*
     * add current point to the heap
     * (add to the bottom, then move up)
     *
     * we don't write point to the heap
     * until its final position is determined
     * (it allow us to reduce number of array access operations)
     */
    j = *n;
    *n = *n+1;
    while(j>0)
    {
        k = (j-1)/2;
        v = a->ptr.p_double[k];
        if( ae_fp_less(v,va) )
        {
            
            /*
             * swap with higher element
             */
            a->ptr.p_double[j] = v;
            b->ptr.p_int[j] = b->ptr.p_int[k];
            j = k;
        }
        else
        {
            
            /*
             * element in its place. terminate.
             */
            break;
        }
    }
    a->ptr.p_double[j] = va;
    b->ptr.p_int[j] = vb;
}


/*************************************************************************
Heap operations: replaces top element with new element
(which is moved down)

PARAMETERS:
    A       -   heap itself, must be at least array[0..N-1]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap
    VA      -   value of the element which replaces top element
    VB      -   value of the tag

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheapreplacetopi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t n,
     double va,
     ae_int_t vb,
     ae_state *_state)
{
    ae_int_t j;
    ae_int_t k1;
    ae_int_t k2;
    double v;
    double v1;
    double v2;


    if( n<1 )
    {
        return;
    }
    
    /*
     * N=1 is a special case
     */
    if( n==1 )
    {
        a->ptr.p_double[0] = va;
        b->ptr.p_int[0] = vb;
        return;
    }
    
    /*
     * move down through heap:
     * * J  -   current element
     * * K1 -   first child (always exists)
     * * K2 -   second child (may not exists)
     *
     * we don't write point to the heap
     * until its final position is determined
     * (it allow us to reduce number of array access operations)
     */
    j = 0;
    k1 = 1;
    k2 = 2;
    while(k1<n)
    {
        if( k2>=n )
        {
            
            /*
             * only one child.
             *
             * swap and terminate (because this child
             * have no siblings due to heap structure)
             */
            v = a->ptr.p_double[k1];
            if( ae_fp_greater(v,va) )
            {
                a->ptr.p_double[j] = v;
                b->ptr.p_int[j] = b->ptr.p_int[k1];
                j = k1;
            }
            break;
        }
        else
        {
            
            /*
             * two childs
             */
            v1 = a->ptr.p_double[k1];
            v2 = a->ptr.p_double[k2];
            if( ae_fp_greater(v1,v2) )
            {
                if( ae_fp_less(va,v1) )
                {
                    a->ptr.p_double[j] = v1;
                    b->ptr.p_int[j] = b->ptr.p_int[k1];
                    j = k1;
                }
                else
                {
                    break;
                }
            }
            else
            {
                if( ae_fp_less(va,v2) )
                {
                    a->ptr.p_double[j] = v2;
                    b->ptr.p_int[j] = b->ptr.p_int[k2];
                    j = k2;
                }
                else
                {
                    break;
                }
            }
            k1 = 2*j+1;
            k2 = 2*j+2;
        }
    }
    a->ptr.p_double[j] = va;
    b->ptr.p_int[j] = vb;
}


/*************************************************************************
Heap operations: pops top element from the heap

PARAMETERS:
    A       -   heap itself, must be at least array[0..N-1]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap, N>=1

On output top element is moved to A[N-1], B[N-1], heap is reordered, N is
decreased by 1.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheappopi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t* n,
     ae_state *_state)
{
    double va;
    ae_int_t vb;


    if( *n<1 )
    {
        return;
    }
    
    /*
     * N=1 is a special case
     */
    if( *n==1 )
    {
        *n = 0;
        return;
    }
    
    /*
     * swap top element and last element,
     * then reorder heap
     */
    va = a->ptr.p_double[*n-1];
    vb = b->ptr.p_int[*n-1];
    a->ptr.p_double[*n-1] = a->ptr.p_double[0];
    b->ptr.p_int[*n-1] = b->ptr.p_int[0];
    *n = *n-1;
    tagheapreplacetopi(a, b, *n, va, vb, _state);
}


/*************************************************************************
Search first element less than T in sorted array.

PARAMETERS:
    A - sorted array by ascending from 0 to N-1
    N - number of elements in array
    T - the desired element

RESULT:
    The very first element's index, which isn't less than T.
In the case when there aren't such elements, returns N.
*************************************************************************/
ae_int_t lowerbound(/* Real    */ ae_vector* a,
     ae_int_t n,
     double t,
     ae_state *_state)
{
    ae_int_t l;
    ae_int_t half;
    ae_int_t first;
    ae_int_t middle;
    ae_int_t result;


    l = n;
    first = 0;
    while(l>0)
    {
        half = l/2;
        middle = first+half;
        if( ae_fp_less(a->ptr.p_double[middle],t) )
        {
            first = middle+1;
            l = l-half-1;
        }
        else
        {
            l = half;
        }
    }
    result = first;
    return result;
}


/*************************************************************************
Search first element more than T in sorted array.

PARAMETERS:
    A - sorted array by ascending from 0 to N-1
    N - number of elements in array
    T - the desired element

    RESULT:
    The very first element's index, which more than T.
In the case when there aren't such elements, returns N.
*************************************************************************/
ae_int_t upperbound(/* Real    */ ae_vector* a,
     ae_int_t n,
     double t,
     ae_state *_state)
{
    ae_int_t l;
    ae_int_t half;
    ae_int_t first;
    ae_int_t middle;
    ae_int_t result;


    l = n;
    first = 0;
    while(l>0)
    {
        half = l/2;
        middle = first+half;
        if( ae_fp_less(t,a->ptr.p_double[middle]) )
        {
            l = half;
        }
        else
        {
            first = middle+1;
            l = l-half-1;
        }
    }
    result = first;
    return result;
}


/*************************************************************************
Internal TagSortFastI: sorts A[I1...I2] (both bounds are included),
applies same permutations to B.

  -- ALGLIB --
     Copyright 06.09.2010 by Bochkanov Sergey
*************************************************************************/
static void tsort_tagsortfastirec(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Integer */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t cntless;
    ae_int_t cnteq;
    ae_int_t cntgreater;
    double tmpr;
    ae_int_t tmpi;
    double v0;
    double v1;
    double v2;
    double vp;


    
    /*
     * Fast exit
     */
    if( i2<=i1 )
    {
        return;
    }
    
    /*
     * Non-recursive sort for small arrays
     */
    if( i2-i1<=16 )
    {
        for(j=i1+1; j<=i2; j++)
        {
            
            /*
             * Search elements [I1..J-1] for place to insert Jth element.
             *
             * This code stops immediately if we can leave A[J] at J-th position
             * (all elements have same value of A[J] larger than any of them)
             */
            tmpr = a->ptr.p_double[j];
            tmpi = j;
            for(k=j-1; k>=i1; k--)
            {
                if( a->ptr.p_double[k]<=tmpr )
                {
                    break;
                }
                tmpi = k;
            }
            k = tmpi;
            
            /*
             * Insert Jth element into Kth position
             */
            if( k!=j )
            {
                tmpr = a->ptr.p_double[j];
                tmpi = b->ptr.p_int[j];
                for(i=j-1; i>=k; i--)
                {
                    a->ptr.p_double[i+1] = a->ptr.p_double[i];
                    b->ptr.p_int[i+1] = b->ptr.p_int[i];
                }
                a->ptr.p_double[k] = tmpr;
                b->ptr.p_int[k] = tmpi;
            }
        }
        return;
    }
    
    /*
     * Quicksort: choose pivot
     * Here we assume that I2-I1>=2
     */
    v0 = a->ptr.p_double[i1];
    v1 = a->ptr.p_double[i1+(i2-i1)/2];
    v2 = a->ptr.p_double[i2];
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    if( v1>v2 )
    {
        tmpr = v2;
        v2 = v1;
        v1 = tmpr;
    }
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    vp = v1;
    
    /*
     * now pass through A/B and:
     * * move elements that are LESS than VP to the left of A/B
     * * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
     * * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
     * * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
     * * move elements from the left of BufA/BufB to the end of A/B
     */
    cntless = 0;
    cnteq = 0;
    cntgreater = 0;
    for(i=i1; i<=i2; i++)
    {
        v0 = a->ptr.p_double[i];
        if( v0<vp )
        {
            
            /*
             * LESS
             */
            k = i1+cntless;
            if( i!=k )
            {
                a->ptr.p_double[k] = v0;
                b->ptr.p_int[k] = b->ptr.p_int[i];
            }
            cntless = cntless+1;
            continue;
        }
        if( v0==vp )
        {
            
            /*
             * EQUAL
             */
            k = i2-cnteq;
            bufa->ptr.p_double[k] = v0;
            bufb->ptr.p_int[k] = b->ptr.p_int[i];
            cnteq = cnteq+1;
            continue;
        }
        
        /*
         * GREATER
         */
        k = i1+cntgreater;
        bufa->ptr.p_double[k] = v0;
        bufb->ptr.p_int[k] = b->ptr.p_int[i];
        cntgreater = cntgreater+1;
    }
    for(i=0; i<=cnteq-1; i++)
    {
        j = i1+cntless+cnteq-1-i;
        k = i2+i-(cnteq-1);
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_int[j] = bufb->ptr.p_int[k];
    }
    for(i=0; i<=cntgreater-1; i++)
    {
        j = i1+cntless+cnteq+i;
        k = i1+i;
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_int[j] = bufb->ptr.p_int[k];
    }
    
    /*
     * Sort left and right parts of the array (ignoring middle part)
     */
    tsort_tagsortfastirec(a, b, bufa, bufb, i1, i1+cntless-1, _state);
    tsort_tagsortfastirec(a, b, bufa, bufb, i1+cntless+cnteq, i2, _state);
}


/*************************************************************************
Internal TagSortFastR: sorts A[I1...I2] (both bounds are included),
applies same permutations to B.

  -- ALGLIB --
     Copyright 06.09.2010 by Bochkanov Sergey
*************************************************************************/
static void tsort_tagsortfastrrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double tmpr;
    double tmpr2;
    ae_int_t tmpi;
    ae_int_t cntless;
    ae_int_t cnteq;
    ae_int_t cntgreater;
    double v0;
    double v1;
    double v2;
    double vp;


    
    /*
     * Fast exit
     */
    if( i2<=i1 )
    {
        return;
    }
    
    /*
     * Non-recursive sort for small arrays
     */
    if( i2-i1<=16 )
    {
        for(j=i1+1; j<=i2; j++)
        {
            
            /*
             * Search elements [I1..J-1] for place to insert Jth element.
             *
             * This code stops immediatly if we can leave A[J] at J-th position
             * (all elements have same value of A[J] larger than any of them)
             */
            tmpr = a->ptr.p_double[j];
            tmpi = j;
            for(k=j-1; k>=i1; k--)
            {
                if( a->ptr.p_double[k]<=tmpr )
                {
                    break;
                }
                tmpi = k;
            }
            k = tmpi;
            
            /*
             * Insert Jth element into Kth position
             */
            if( k!=j )
            {
                tmpr = a->ptr.p_double[j];
                tmpr2 = b->ptr.p_double[j];
                for(i=j-1; i>=k; i--)
                {
                    a->ptr.p_double[i+1] = a->ptr.p_double[i];
                    b->ptr.p_double[i+1] = b->ptr.p_double[i];
                }
                a->ptr.p_double[k] = tmpr;
                b->ptr.p_double[k] = tmpr2;
            }
        }
        return;
    }
    
    /*
     * Quicksort: choose pivot
     * Here we assume that I2-I1>=16
     */
    v0 = a->ptr.p_double[i1];
    v1 = a->ptr.p_double[i1+(i2-i1)/2];
    v2 = a->ptr.p_double[i2];
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    if( v1>v2 )
    {
        tmpr = v2;
        v2 = v1;
        v1 = tmpr;
    }
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    vp = v1;
    
    /*
     * now pass through A/B and:
     * * move elements that are LESS than VP to the left of A/B
     * * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
     * * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
     * * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
     * * move elements from the left of BufA/BufB to the end of A/B
     */
    cntless = 0;
    cnteq = 0;
    cntgreater = 0;
    for(i=i1; i<=i2; i++)
    {
        v0 = a->ptr.p_double[i];
        if( v0<vp )
        {
            
            /*
             * LESS
             */
            k = i1+cntless;
            if( i!=k )
            {
                a->ptr.p_double[k] = v0;
                b->ptr.p_double[k] = b->ptr.p_double[i];
            }
            cntless = cntless+1;
            continue;
        }
        if( v0==vp )
        {
            
            /*
             * EQUAL
             */
            k = i2-cnteq;
            bufa->ptr.p_double[k] = v0;
            bufb->ptr.p_double[k] = b->ptr.p_double[i];
            cnteq = cnteq+1;
            continue;
        }
        
        /*
         * GREATER
         */
        k = i1+cntgreater;
        bufa->ptr.p_double[k] = v0;
        bufb->ptr.p_double[k] = b->ptr.p_double[i];
        cntgreater = cntgreater+1;
    }
    for(i=0; i<=cnteq-1; i++)
    {
        j = i1+cntless+cnteq-1-i;
        k = i2+i-(cnteq-1);
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_double[j] = bufb->ptr.p_double[k];
    }
    for(i=0; i<=cntgreater-1; i++)
    {
        j = i1+cntless+cnteq+i;
        k = i1+i;
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_double[j] = bufb->ptr.p_double[k];
    }
    
    /*
     * Sort left and right parts of the array (ignoring middle part)
     */
    tsort_tagsortfastrrec(a, b, bufa, bufb, i1, i1+cntless-1, _state);
    tsort_tagsortfastrrec(a, b, bufa, bufb, i1+cntless+cnteq, i2, _state);
}


/*************************************************************************
Internal TagSortFastI: sorts A[I1...I2] (both bounds are included),
applies same permutations to B.

  -- ALGLIB --
     Copyright 06.09.2010 by Bochkanov Sergey
*************************************************************************/
static void tsort_tagsortfastrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* bufa,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t cntless;
    ae_int_t cnteq;
    ae_int_t cntgreater;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double tmpr;
    ae_int_t tmpi;
    double v0;
    double v1;
    double v2;
    double vp;


    
    /*
     * Fast exit
     */
    if( i2<=i1 )
    {
        return;
    }
    
    /*
     * Non-recursive sort for small arrays
     */
    if( i2-i1<=16 )
    {
        for(j=i1+1; j<=i2; j++)
        {
            
            /*
             * Search elements [I1..J-1] for place to insert Jth element.
             *
             * This code stops immediatly if we can leave A[J] at J-th position
             * (all elements have same value of A[J] larger than any of them)
             */
            tmpr = a->ptr.p_double[j];
            tmpi = j;
            for(k=j-1; k>=i1; k--)
            {
                if( a->ptr.p_double[k]<=tmpr )
                {
                    break;
                }
                tmpi = k;
            }
            k = tmpi;
            
            /*
             * Insert Jth element into Kth position
             */
            if( k!=j )
            {
                tmpr = a->ptr.p_double[j];
                for(i=j-1; i>=k; i--)
                {
                    a->ptr.p_double[i+1] = a->ptr.p_double[i];
                }
                a->ptr.p_double[k] = tmpr;
            }
        }
        return;
    }
    
    /*
     * Quicksort: choose pivot
     * Here we assume that I2-I1>=16
     */
    v0 = a->ptr.p_double[i1];
    v1 = a->ptr.p_double[i1+(i2-i1)/2];
    v2 = a->ptr.p_double[i2];
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    if( v1>v2 )
    {
        tmpr = v2;
        v2 = v1;
        v1 = tmpr;
    }
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    vp = v1;
    
    /*
     * now pass through A/B and:
     * * move elements that are LESS than VP to the left of A/B
     * * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
     * * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
     * * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
     * * move elements from the left of BufA/BufB to the end of A/B
     */
    cntless = 0;
    cnteq = 0;
    cntgreater = 0;
    for(i=i1; i<=i2; i++)
    {
        v0 = a->ptr.p_double[i];
        if( v0<vp )
        {
            
            /*
             * LESS
             */
            k = i1+cntless;
            if( i!=k )
            {
                a->ptr.p_double[k] = v0;
            }
            cntless = cntless+1;
            continue;
        }
        if( v0==vp )
        {
            
            /*
             * EQUAL
             */
            k = i2-cnteq;
            bufa->ptr.p_double[k] = v0;
            cnteq = cnteq+1;
            continue;
        }
        
        /*
         * GREATER
         */
        k = i1+cntgreater;
        bufa->ptr.p_double[k] = v0;
        cntgreater = cntgreater+1;
    }
    for(i=0; i<=cnteq-1; i++)
    {
        j = i1+cntless+cnteq-1-i;
        k = i2+i-(cnteq-1);
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
    }
    for(i=0; i<=cntgreater-1; i++)
    {
        j = i1+cntless+cnteq+i;
        k = i1+i;
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
    }
    
    /*
     * Sort left and right parts of the array (ignoring middle part)
     */
    tsort_tagsortfastrec(a, bufa, i1, i1+cntless-1, _state);
    tsort_tagsortfastrec(a, bufa, i1+cntless+cnteq, i2, _state);
}




/*************************************************************************
Internal ranking subroutine.

INPUT PARAMETERS:
    X       -   array to rank
    N       -   array size
    IsCentered- whether ranks are centered or not:
                * True      -   ranks are centered in such way that  their
                                sum is zero
                * False     -   ranks are not centered
    Buf     -   temporary buffers
    
NOTE: when IsCentered is True and all X[] are equal, this  function  fills
      X by zeros (exact zeros are used, not sum which is only approximately
      equal to zero).
*************************************************************************/
void rankx(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_bool iscentered,
     apbuffers* buf,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double tmp;
    double voffs;


    
    /*
     * Prepare
     */
    if( n<1 )
    {
        return;
    }
    if( n==1 )
    {
        x->ptr.p_double[0] = (double)(0);
        return;
    }
    if( buf->ra1.cnt<n )
    {
        ae_vector_set_length(&buf->ra1, n, _state);
    }
    if( buf->ia1.cnt<n )
    {
        ae_vector_set_length(&buf->ia1, n, _state);
    }
    for(i=0; i<=n-1; i++)
    {
        buf->ra1.ptr.p_double[i] = x->ptr.p_double[i];
        buf->ia1.ptr.p_int[i] = i;
    }
    tagsortfasti(&buf->ra1, &buf->ia1, &buf->ra2, &buf->ia2, n, _state);
    
    /*
     * Special test for all values being equal
     */
    if( ae_fp_eq(buf->ra1.ptr.p_double[0],buf->ra1.ptr.p_double[n-1]) )
    {
        if( iscentered )
        {
            tmp = 0.0;
        }
        else
        {
            tmp = (double)(n-1)/(double)2;
        }
        for(i=0; i<=n-1; i++)
        {
            x->ptr.p_double[i] = tmp;
        }
        return;
    }
    
    /*
     * compute tied ranks
     */
    i = 0;
    while(i<=n-1)
    {
        j = i+1;
        while(j<=n-1)
        {
            if( ae_fp_neq(buf->ra1.ptr.p_double[j],buf->ra1.ptr.p_double[i]) )
            {
                break;
            }
            j = j+1;
        }
        for(k=i; k<=j-1; k++)
        {
            buf->ra1.ptr.p_double[k] = (double)(i+j-1)/(double)2;
        }
        i = j;
    }
    
    /*
     * back to x
     */
    if( iscentered )
    {
        voffs = (double)(n-1)/(double)2;
    }
    else
    {
        voffs = 0.0;
    }
    for(i=0; i<=n-1; i++)
    {
        x->ptr.p_double[buf->ia1.ptr.p_int[i]] = buf->ra1.ptr.p_double[i]-voffs;
    }
}




/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixrank1f(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     /* Complex */ ae_vector* u,
     ae_int_t iu,
     /* Complex */ ae_vector* v,
     ae_int_t iv,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixrank1f(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     /* Real    */ ae_vector* u,
     ae_int_t iu,
     /* Real    */ ae_vector* v,
     ae_int_t iv,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixmvf(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t opa,
     /* Complex */ ae_vector* x,
     ae_int_t ix,
     /* Complex */ ae_vector* y,
     ae_int_t iy,
     ae_state *_state)
{
    ae_bool result;


    result = ae_false;
    return result;
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixmvf(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t opa,
     /* Real    */ ae_vector* x,
     ae_int_t ix,
     /* Real    */ ae_vector* y,
     ae_int_t iy,
     ae_state *_state)
{
    ae_bool result;


    result = ae_false;
    return result;
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixrighttrsmf(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Complex */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixlefttrsmf(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Complex */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixrighttrsmf(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Real    */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixlefttrsmf(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Real    */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixherkf(ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     double beta,
     /* Complex */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_bool isupper,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixherkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixsyrkf(ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_bool isupper,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixsyrkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixgemmf(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     /* Real    */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     ae_int_t optypeb,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixgemmf(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     ae_complex alpha,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     /* Complex */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     ae_int_t optypeb,
     ae_complex beta,
     /* Complex */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}


/*************************************************************************
CMatrixGEMM kernel, basecase code for CMatrixGEMM.

This subroutine calculates C = alpha*op1(A)*op2(B) +beta*C where:
* C is MxN general matrix
* op1(A) is MxK matrix
* op2(B) is KxN matrix
* "op" may be identity transformation, transposition, conjugate transposition

Additional info:
* multiplication result replaces C. If Beta=0, C elements are not used in
  calculations (not multiplied by zero - just not referenced)
* if Alpha=0, A is not used (not multiplied by zero - just not referenced)
* if both Beta and Alpha are zero, C is filled by zeros.

IMPORTANT:

This function does NOT preallocate output matrix C, it MUST be preallocated
by caller prior to calling this function. In case C does not have  enough
space to store result, exception will be generated.

INPUT PARAMETERS
    M       -   matrix size, M>0
    N       -   matrix size, N>0
    K       -   matrix size, K>0
    Alpha   -   coefficient
    A       -   matrix
    IA      -   submatrix offset
    JA      -   submatrix offset
    OpTypeA -   transformation type:
                * 0 - no transformation
                * 1 - transposition
                * 2 - conjugate transposition
    B       -   matrix
    IB      -   submatrix offset
    JB      -   submatrix offset
    OpTypeB -   transformation type:
                * 0 - no transformation
                * 1 - transposition
                * 2 - conjugate transposition
    Beta    -   coefficient
    C       -   PREALLOCATED output matrix
    IC      -   submatrix offset
    JC      -   submatrix offset

  -- ALGLIB routine --
     27.03.2013
     Bochkanov Sergey
*************************************************************************/
void cmatrixgemmk(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     ae_complex alpha,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     /* Complex */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     ae_int_t optypeb,
     ae_complex beta,
     /* Complex */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_complex v;
    ae_complex v00;
    ae_complex v01;
    ae_complex v10;
    ae_complex v11;
    double v00x;
    double v00y;
    double v01x;
    double v01y;
    double v10x;
    double v10y;
    double v11x;
    double v11y;
    double a0x;
    double a0y;
    double a1x;
    double a1y;
    double b0x;
    double b0y;
    double b1x;
    double b1y;
    ae_int_t idxa0;
    ae_int_t idxa1;
    ae_int_t idxb0;
    ae_int_t idxb1;
    ae_int_t i0;
    ae_int_t i1;
    ae_int_t ik;
    ae_int_t j0;
    ae_int_t j1;
    ae_int_t jk;
    ae_int_t t;
    ae_int_t offsa;
    ae_int_t offsb;


    
    /*
     * if matrix size is zero
     */
    if( m==0||n==0 )
    {
        return;
    }
    
    /*
     * Try optimized code
     */
    if( cmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc, _state) )
    {
        return;
    }
    
    /*
     * if K=0, then C=Beta*C
     */
    if( k==0 )
    {
        if( ae_c_neq_d(beta,(double)(1)) )
        {
            if( ae_c_neq_d(beta,(double)(0)) )
            {
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c->ptr.pp_complex[ic+i][jc+j] = ae_c_mul(beta,c->ptr.pp_complex[ic+i][jc+j]);
                    }
                }
            }
            else
            {
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c->ptr.pp_complex[ic+i][jc+j] = ae_complex_from_i(0);
                    }
                }
            }
        }
        return;
    }
    
    /*
     * This phase is not really necessary, but compiler complains
     * about "possibly uninitialized variables"
     */
    a0x = (double)(0);
    a0y = (double)(0);
    a1x = (double)(0);
    a1y = (double)(0);
    b0x = (double)(0);
    b0y = (double)(0);
    b1x = (double)(0);
    b1y = (double)(0);
    
    /*
     * General case
     */
    i = 0;
    while(i<m)
    {
        j = 0;
        while(j<n)
        {
            
            /*
             * Choose between specialized 4x4 code and general code
             */
            if( i+2<=m&&j+2<=n )
            {
                
                /*
                 * Specialized 4x4 code for [I..I+3]x[J..J+3] submatrix of C.
                 *
                 * This submatrix is calculated as sum of K rank-1 products,
                 * with operands cached in local variables in order to speed
                 * up operations with arrays.
                 */
                v00x = 0.0;
                v00y = 0.0;
                v01x = 0.0;
                v01y = 0.0;
                v10x = 0.0;
                v10y = 0.0;
                v11x = 0.0;
                v11y = 0.0;
                if( optypea==0 )
                {
                    idxa0 = ia+i+0;
                    idxa1 = ia+i+1;
                    offsa = ja;
                }
                else
                {
                    idxa0 = ja+i+0;
                    idxa1 = ja+i+1;
                    offsa = ia;
                }
                if( optypeb==0 )
                {
                    idxb0 = jb+j+0;
                    idxb1 = jb+j+1;
                    offsb = ib;
                }
                else
                {
                    idxb0 = ib+j+0;
                    idxb1 = ib+j+1;
                    offsb = jb;
                }
                for(t=0; t<=k-1; t++)
                {
                    if( optypea==0 )
                    {
                        a0x = a->ptr.pp_complex[idxa0][offsa].x;
                        a0y = a->ptr.pp_complex[idxa0][offsa].y;
                        a1x = a->ptr.pp_complex[idxa1][offsa].x;
                        a1y = a->ptr.pp_complex[idxa1][offsa].y;
                    }
                    if( optypea==1 )
                    {
                        a0x = a->ptr.pp_complex[offsa][idxa0].x;
                        a0y = a->ptr.pp_complex[offsa][idxa0].y;
                        a1x = a->ptr.pp_complex[offsa][idxa1].x;
                        a1y = a->ptr.pp_complex[offsa][idxa1].y;
                    }
                    if( optypea==2 )
                    {
                        a0x = a->ptr.pp_complex[offsa][idxa0].x;
                        a0y = -a->ptr.pp_complex[offsa][idxa0].y;
                        a1x = a->ptr.pp_complex[offsa][idxa1].x;
                        a1y = -a->ptr.pp_complex[offsa][idxa1].y;
                    }
                    if( optypeb==0 )
                    {
                        b0x = b->ptr.pp_complex[offsb][idxb0].x;
                        b0y = b->ptr.pp_complex[offsb][idxb0].y;
                        b1x = b->ptr.pp_complex[offsb][idxb1].x;
                        b1y = b->ptr.pp_complex[offsb][idxb1].y;
                    }
                    if( optypeb==1 )
                    {
                        b0x = b->ptr.pp_complex[idxb0][offsb].x;
                        b0y = b->ptr.pp_complex[idxb0][offsb].y;
                        b1x = b->ptr.pp_complex[idxb1][offsb].x;
                        b1y = b->ptr.pp_complex[idxb1][offsb].y;
                    }
                    if( optypeb==2 )
                    {
                        b0x = b->ptr.pp_complex[idxb0][offsb].x;
                        b0y = -b->ptr.pp_complex[idxb0][offsb].y;
                        b1x = b->ptr.pp_complex[idxb1][offsb].x;
                        b1y = -b->ptr.pp_complex[idxb1][offsb].y;
                    }
                    v00x = v00x+a0x*b0x-a0y*b0y;
                    v00y = v00y+a0x*b0y+a0y*b0x;
                    v01x = v01x+a0x*b1x-a0y*b1y;
                    v01y = v01y+a0x*b1y+a0y*b1x;
                    v10x = v10x+a1x*b0x-a1y*b0y;
                    v10y = v10y+a1x*b0y+a1y*b0x;
                    v11x = v11x+a1x*b1x-a1y*b1y;
                    v11y = v11y+a1x*b1y+a1y*b1x;
                    offsa = offsa+1;
                    offsb = offsb+1;
                }
                v00.x = v00x;
                v00.y = v00y;
                v10.x = v10x;
                v10.y = v10y;
                v01.x = v01x;
                v01.y = v01y;
                v11.x = v11x;
                v11.y = v11y;
                if( ae_c_eq_d(beta,(double)(0)) )
                {
                    c->ptr.pp_complex[ic+i+0][jc+j+0] = ae_c_mul(alpha,v00);
                    c->ptr.pp_complex[ic+i+0][jc+j+1] = ae_c_mul(alpha,v01);
                    c->ptr.pp_complex[ic+i+1][jc+j+0] = ae_c_mul(alpha,v10);
                    c->ptr.pp_complex[ic+i+1][jc+j+1] = ae_c_mul(alpha,v11);
                }
                else
                {
                    c->ptr.pp_complex[ic+i+0][jc+j+0] = ae_c_add(ae_c_mul(beta,c->ptr.pp_complex[ic+i+0][jc+j+0]),ae_c_mul(alpha,v00));
                    c->ptr.pp_complex[ic+i+0][jc+j+1] = ae_c_add(ae_c_mul(beta,c->ptr.pp_complex[ic+i+0][jc+j+1]),ae_c_mul(alpha,v01));
                    c->ptr.pp_complex[ic+i+1][jc+j+0] = ae_c_add(ae_c_mul(beta,c->ptr.pp_complex[ic+i+1][jc+j+0]),ae_c_mul(alpha,v10));
                    c->ptr.pp_complex[ic+i+1][jc+j+1] = ae_c_add(ae_c_mul(beta,c->ptr.pp_complex[ic+i+1][jc+j+1]),ae_c_mul(alpha,v11));
                }
            }
            else
            {
                
                /*
                 * Determine submatrix [I0..I1]x[J0..J1] to process
                 */
                i0 = i;
                i1 = ae_minint(i+1, m-1, _state);
                j0 = j;
                j1 = ae_minint(j+1, n-1, _state);
                
                /*
                 * Process submatrix
                 */
                for(ik=i0; ik<=i1; ik++)
                {
                    for(jk=j0; jk<=j1; jk++)
                    {
                        if( k==0||ae_c_eq_d(alpha,(double)(0)) )
                        {
                            v = ae_complex_from_i(0);
                        }
                        else
                        {
                            v = ae_complex_from_d(0.0);
                            if( optypea==0&&optypeb==0 )
                            {
                                v = ae_v_cdotproduct(&a->ptr.pp_complex[ia+ik][ja], 1, "N", &b->ptr.pp_complex[ib][jb+jk], b->stride, "N", ae_v_len(ja,ja+k-1));
                            }
                            if( optypea==0&&optypeb==1 )
                            {
                                v = ae_v_cdotproduct(&a->ptr.pp_complex[ia+ik][ja], 1, "N", &b->ptr.pp_complex[ib+jk][jb], 1, "N", ae_v_len(ja,ja+k-1));
                            }
                            if( optypea==0&&optypeb==2 )
                            {
                                v = ae_v_cdotproduct(&a->ptr.pp_complex[ia+ik][ja], 1, "N", &b->ptr.pp_complex[ib+jk][jb], 1, "Conj", ae_v_len(ja,ja+k-1));
                            }
                            if( optypea==1&&optypeb==0 )
                            {
                                v = ae_v_cdotproduct(&a->ptr.pp_complex[ia][ja+ik], a->stride, "N", &b->ptr.pp_complex[ib][jb+jk], b->stride, "N", ae_v_len(ia,ia+k-1));
                            }
                            if( optypea==1&&optypeb==1 )
                            {
                                v = ae_v_cdotproduct(&a->ptr.pp_complex[ia][ja+ik], a->stride, "N", &b->ptr.pp_complex[ib+jk][jb], 1, "N", ae_v_len(ia,ia+k-1));
                            }
                            if( optypea==1&&optypeb==2 )
                            {
                                v = ae_v_cdotproduct(&a->ptr.pp_complex[ia][ja+ik], a->stride, "N", &b->ptr.pp_complex[ib+jk][jb], 1, "Conj", ae_v_len(ia,ia+k-1));
                            }
                            if( optypea==2&&optypeb==0 )
                            {
                                v = ae_v_cdotproduct(&a->ptr.pp_complex[ia][ja+ik], a->stride, "Conj", &b->ptr.pp_complex[ib][jb+jk], b->stride, "N", ae_v_len(ia,ia+k-1));
                            }
                            if( optypea==2&&optypeb==1 )
                            {
                                v = ae_v_cdotproduct(&a->ptr.pp_complex[ia][ja+ik], a->stride, "Conj", &b->ptr.pp_complex[ib+jk][jb], 1, "N", ae_v_len(ia,ia+k-1));
                            }
                            if( optypea==2&&optypeb==2 )
                            {
                                v = ae_v_cdotproduct(&a->ptr.pp_complex[ia][ja+ik], a->stride, "Conj", &b->ptr.pp_complex[ib+jk][jb], 1, "Conj", ae_v_len(ia,ia+k-1));
                            }
                        }
                        if( ae_c_eq_d(beta,(double)(0)) )
                        {
                            c->ptr.pp_complex[ic+ik][jc+jk] = ae_c_mul(alpha,v);
                        }
                        else
                        {
                            c->ptr.pp_complex[ic+ik][jc+jk] = ae_c_add(ae_c_mul(beta,c->ptr.pp_complex[ic+ik][jc+jk]),ae_c_mul(alpha,v));
                        }
                    }
                }
            }
            j = j+2;
        }
        i = i+2;
    }
}


/*************************************************************************
RMatrixGEMM kernel, basecase code for RMatrixGEMM.

This subroutine calculates C = alpha*op1(A)*op2(B) +beta*C where:
* C is MxN general matrix
* op1(A) is MxK matrix
* op2(B) is KxN matrix
* "op" may be identity transformation, transposition

Additional info:
* multiplication result replaces C. If Beta=0, C elements are not used in
  calculations (not multiplied by zero - just not referenced)
* if Alpha=0, A is not used (not multiplied by zero - just not referenced)
* if both Beta and Alpha are zero, C is filled by zeros.

IMPORTANT:

This function does NOT preallocate output matrix C, it MUST be preallocated
by caller prior to calling this function. In case C does not have  enough
space to store result, exception will be generated.

INPUT PARAMETERS
    M       -   matrix size, M>0
    N       -   matrix size, N>0
    K       -   matrix size, K>0
    Alpha   -   coefficient
    A       -   matrix
    IA      -   submatrix offset
    JA      -   submatrix offset
    OpTypeA -   transformation type:
                * 0 - no transformation
                * 1 - transposition
    B       -   matrix
    IB      -   submatrix offset
    JB      -   submatrix offset
    OpTypeB -   transformation type:
                * 0 - no transformation
                * 1 - transposition
    Beta    -   coefficient
    C       -   PREALLOCATED output matrix
    IC      -   submatrix offset
    JC      -   submatrix offset

  -- ALGLIB routine --
     27.03.2013
     Bochkanov Sergey
*************************************************************************/
void rmatrixgemmk(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     /* Real    */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     ae_int_t optypeb,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    
    /*
     * if matrix size is zero
     */
    if( m==0||n==0 )
    {
        return;
    }
    
    /*
     * Try optimized code
     */
    if( rmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc, _state) )
    {
        return;
    }
    
    /*
     * if K=0, then C=Beta*C
     */
    if( k==0||ae_fp_eq(alpha,(double)(0)) )
    {
        if( ae_fp_neq(beta,(double)(1)) )
        {
            if( ae_fp_neq(beta,(double)(0)) )
            {
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c->ptr.pp_double[ic+i][jc+j] = beta*c->ptr.pp_double[ic+i][jc+j];
                    }
                }
            }
            else
            {
                for(i=0; i<=m-1; i++)
                {
                    for(j=0; j<=n-1; j++)
                    {
                        c->ptr.pp_double[ic+i][jc+j] = (double)(0);
                    }
                }
            }
        }
        return;
    }
    
    /*
     * Call specialized code.
     *
     * NOTE: specialized code was moved to separate function because of strange
     *       issues with instructions cache on some systems; Having too long
     *       functions significantly slows down internal loop of the algorithm.
     */
    if( optypea==0&&optypeb==0 )
    {
        rmatrixgemmk44v00(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc, _state);
    }
    if( optypea==0&&optypeb!=0 )
    {
        rmatrixgemmk44v01(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc, _state);
    }
    if( optypea!=0&&optypeb==0 )
    {
        rmatrixgemmk44v10(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc, _state);
    }
    if( optypea!=0&&optypeb!=0 )
    {
        rmatrixgemmk44v11(m, n, k, alpha, a, ia, ja, b, ib, jb, beta, c, ic, jc, _state);
    }
}


/*************************************************************************
RMatrixGEMM kernel, basecase code for RMatrixGEMM, specialized for sitation
with OpTypeA=0 and OpTypeB=0.

Additional info:
* this function requires that Alpha<>0 (assertion is thrown otherwise)

INPUT PARAMETERS
    M       -   matrix size, M>0
    N       -   matrix size, N>0
    K       -   matrix size, K>0
    Alpha   -   coefficient
    A       -   matrix
    IA      -   submatrix offset
    JA      -   submatrix offset
    B       -   matrix
    IB      -   submatrix offset
    JB      -   submatrix offset
    Beta    -   coefficient
    C       -   PREALLOCATED output matrix
    IC      -   submatrix offset
    JC      -   submatrix offset

  -- ALGLIB routine --
     27.03.2013
     Bochkanov Sergey
*************************************************************************/
void rmatrixgemmk44v00(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     /* Real    */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double v;
    double v00;
    double v01;
    double v02;
    double v03;
    double v10;
    double v11;
    double v12;
    double v13;
    double v20;
    double v21;
    double v22;
    double v23;
    double v30;
    double v31;
    double v32;
    double v33;
    double a0;
    double a1;
    double a2;
    double a3;
    double b0;
    double b1;
    double b2;
    double b3;
    ae_int_t idxa0;
    ae_int_t idxa1;
    ae_int_t idxa2;
    ae_int_t idxa3;
    ae_int_t idxb0;
    ae_int_t idxb1;
    ae_int_t idxb2;
    ae_int_t idxb3;
    ae_int_t i0;
    ae_int_t i1;
    ae_int_t ik;
    ae_int_t j0;
    ae_int_t j1;
    ae_int_t jk;
    ae_int_t t;
    ae_int_t offsa;
    ae_int_t offsb;


    ae_assert(ae_fp_neq(alpha,(double)(0)), "RMatrixGEMMK44V00: internal error (Alpha=0)", _state);
    
    /*
     * if matrix size is zero
     */
    if( m==0||n==0 )
    {
        return;
    }
    
    /*
     * A*B
     */
    i = 0;
    while(i<m)
    {
        j = 0;
        while(j<n)
        {
            
            /*
             * Choose between specialized 4x4 code and general code
             */
            if( i+4<=m&&j+4<=n )
            {
                
                /*
                 * Specialized 4x4 code for [I..I+3]x[J..J+3] submatrix of C.
                 *
                 * This submatrix is calculated as sum of K rank-1 products,
                 * with operands cached in local variables in order to speed
                 * up operations with arrays.
                 */
                idxa0 = ia+i+0;
                idxa1 = ia+i+1;
                idxa2 = ia+i+2;
                idxa3 = ia+i+3;
                offsa = ja;
                idxb0 = jb+j+0;
                idxb1 = jb+j+1;
                idxb2 = jb+j+2;
                idxb3 = jb+j+3;
                offsb = ib;
                v00 = 0.0;
                v01 = 0.0;
                v02 = 0.0;
                v03 = 0.0;
                v10 = 0.0;
                v11 = 0.0;
                v12 = 0.0;
                v13 = 0.0;
                v20 = 0.0;
                v21 = 0.0;
                v22 = 0.0;
                v23 = 0.0;
                v30 = 0.0;
                v31 = 0.0;
                v32 = 0.0;
                v33 = 0.0;
                
                /*
                 * Different variants of internal loop
                 */
                for(t=0; t<=k-1; t++)
                {
                    a0 = a->ptr.pp_double[idxa0][offsa];
                    a1 = a->ptr.pp_double[idxa1][offsa];
                    b0 = b->ptr.pp_double[offsb][idxb0];
                    b1 = b->ptr.pp_double[offsb][idxb1];
                    v00 = v00+a0*b0;
                    v01 = v01+a0*b1;
                    v10 = v10+a1*b0;
                    v11 = v11+a1*b1;
                    a2 = a->ptr.pp_double[idxa2][offsa];
                    a3 = a->ptr.pp_double[idxa3][offsa];
                    v20 = v20+a2*b0;
                    v21 = v21+a2*b1;
                    v30 = v30+a3*b0;
                    v31 = v31+a3*b1;
                    b2 = b->ptr.pp_double[offsb][idxb2];
                    b3 = b->ptr.pp_double[offsb][idxb3];
                    v22 = v22+a2*b2;
                    v23 = v23+a2*b3;
                    v32 = v32+a3*b2;
                    v33 = v33+a3*b3;
                    v02 = v02+a0*b2;
                    v03 = v03+a0*b3;
                    v12 = v12+a1*b2;
                    v13 = v13+a1*b3;
                    offsa = offsa+1;
                    offsb = offsb+1;
                }
                if( ae_fp_eq(beta,(double)(0)) )
                {
                    c->ptr.pp_double[ic+i+0][jc+j+0] = alpha*v00;
                    c->ptr.pp_double[ic+i+0][jc+j+1] = alpha*v01;
                    c->ptr.pp_double[ic+i+0][jc+j+2] = alpha*v02;
                    c->ptr.pp_double[ic+i+0][jc+j+3] = alpha*v03;
                    c->ptr.pp_double[ic+i+1][jc+j+0] = alpha*v10;
                    c->ptr.pp_double[ic+i+1][jc+j+1] = alpha*v11;
                    c->ptr.pp_double[ic+i+1][jc+j+2] = alpha*v12;
                    c->ptr.pp_double[ic+i+1][jc+j+3] = alpha*v13;
                    c->ptr.pp_double[ic+i+2][jc+j+0] = alpha*v20;
                    c->ptr.pp_double[ic+i+2][jc+j+1] = alpha*v21;
                    c->ptr.pp_double[ic+i+2][jc+j+2] = alpha*v22;
                    c->ptr.pp_double[ic+i+2][jc+j+3] = alpha*v23;
                    c->ptr.pp_double[ic+i+3][jc+j+0] = alpha*v30;
                    c->ptr.pp_double[ic+i+3][jc+j+1] = alpha*v31;
                    c->ptr.pp_double[ic+i+3][jc+j+2] = alpha*v32;
                    c->ptr.pp_double[ic+i+3][jc+j+3] = alpha*v33;
                }
                else
                {
                    c->ptr.pp_double[ic+i+0][jc+j+0] = beta*c->ptr.pp_double[ic+i+0][jc+j+0]+alpha*v00;
                    c->ptr.pp_double[ic+i+0][jc+j+1] = beta*c->ptr.pp_double[ic+i+0][jc+j+1]+alpha*v01;
                    c->ptr.pp_double[ic+i+0][jc+j+2] = beta*c->ptr.pp_double[ic+i+0][jc+j+2]+alpha*v02;
                    c->ptr.pp_double[ic+i+0][jc+j+3] = beta*c->ptr.pp_double[ic+i+0][jc+j+3]+alpha*v03;
                    c->ptr.pp_double[ic+i+1][jc+j+0] = beta*c->ptr.pp_double[ic+i+1][jc+j+0]+alpha*v10;
                    c->ptr.pp_double[ic+i+1][jc+j+1] = beta*c->ptr.pp_double[ic+i+1][jc+j+1]+alpha*v11;
                    c->ptr.pp_double[ic+i+1][jc+j+2] = beta*c->ptr.pp_double[ic+i+1][jc+j+2]+alpha*v12;
                    c->ptr.pp_double[ic+i+1][jc+j+3] = beta*c->ptr.pp_double[ic+i+1][jc+j+3]+alpha*v13;
                    c->ptr.pp_double[ic+i+2][jc+j+0] = beta*c->ptr.pp_double[ic+i+2][jc+j+0]+alpha*v20;
                    c->ptr.pp_double[ic+i+2][jc+j+1] = beta*c->ptr.pp_double[ic+i+2][jc+j+1]+alpha*v21;
                    c->ptr.pp_double[ic+i+2][jc+j+2] = beta*c->ptr.pp_double[ic+i+2][jc+j+2]+alpha*v22;
                    c->ptr.pp_double[ic+i+2][jc+j+3] = beta*c->ptr.pp_double[ic+i+2][jc+j+3]+alpha*v23;
                    c->ptr.pp_double[ic+i+3][jc+j+0] = beta*c->ptr.pp_double[ic+i+3][jc+j+0]+alpha*v30;
                    c->ptr.pp_double[ic+i+3][jc+j+1] = beta*c->ptr.pp_double[ic+i+3][jc+j+1]+alpha*v31;
                    c->ptr.pp_double[ic+i+3][jc+j+2] = beta*c->ptr.pp_double[ic+i+3][jc+j+2]+alpha*v32;
                    c->ptr.pp_double[ic+i+3][jc+j+3] = beta*c->ptr.pp_double[ic+i+3][jc+j+3]+alpha*v33;
                }
            }
            else
            {
                
                /*
                 * Determine submatrix [I0..I1]x[J0..J1] to process
                 */
                i0 = i;
                i1 = ae_minint(i+3, m-1, _state);
                j0 = j;
                j1 = ae_minint(j+3, n-1, _state);
                
                /*
                 * Process submatrix
                 */
                for(ik=i0; ik<=i1; ik++)
                {
                    for(jk=j0; jk<=j1; jk++)
                    {
                        if( k==0||ae_fp_eq(alpha,(double)(0)) )
                        {
                            v = (double)(0);
                        }
                        else
                        {
                            v = ae_v_dotproduct(&a->ptr.pp_double[ia+ik][ja], 1, &b->ptr.pp_double[ib][jb+jk], b->stride, ae_v_len(ja,ja+k-1));
                        }
                        if( ae_fp_eq(beta,(double)(0)) )
                        {
                            c->ptr.pp_double[ic+ik][jc+jk] = alpha*v;
                        }
                        else
                        {
                            c->ptr.pp_double[ic+ik][jc+jk] = beta*c->ptr.pp_double[ic+ik][jc+jk]+alpha*v;
                        }
                    }
                }
            }
            j = j+4;
        }
        i = i+4;
    }
}


/*************************************************************************
RMatrixGEMM kernel, basecase code for RMatrixGEMM, specialized for sitation
with OpTypeA=0 and OpTypeB=1.

Additional info:
* this function requires that Alpha<>0 (assertion is thrown otherwise)

INPUT PARAMETERS
    M       -   matrix size, M>0
    N       -   matrix size, N>0
    K       -   matrix size, K>0
    Alpha   -   coefficient
    A       -   matrix
    IA      -   submatrix offset
    JA      -   submatrix offset
    B       -   matrix
    IB      -   submatrix offset
    JB      -   submatrix offset
    Beta    -   coefficient
    C       -   PREALLOCATED output matrix
    IC      -   submatrix offset
    JC      -   submatrix offset

  -- ALGLIB routine --
     27.03.2013
     Bochkanov Sergey
*************************************************************************/
void rmatrixgemmk44v01(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     /* Real    */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double v;
    double v00;
    double v01;
    double v02;
    double v03;
    double v10;
    double v11;
    double v12;
    double v13;
    double v20;
    double v21;
    double v22;
    double v23;
    double v30;
    double v31;
    double v32;
    double v33;
    double a0;
    double a1;
    double a2;
    double a3;
    double b0;
    double b1;
    double b2;
    double b3;
    ae_int_t idxa0;
    ae_int_t idxa1;
    ae_int_t idxa2;
    ae_int_t idxa3;
    ae_int_t idxb0;
    ae_int_t idxb1;
    ae_int_t idxb2;
    ae_int_t idxb3;
    ae_int_t i0;
    ae_int_t i1;
    ae_int_t ik;
    ae_int_t j0;
    ae_int_t j1;
    ae_int_t jk;
    ae_int_t t;
    ae_int_t offsa;
    ae_int_t offsb;


    ae_assert(ae_fp_neq(alpha,(double)(0)), "RMatrixGEMMK44V00: internal error (Alpha=0)", _state);
    
    /*
     * if matrix size is zero
     */
    if( m==0||n==0 )
    {
        return;
    }
    
    /*
     * A*B'
     */
    i = 0;
    while(i<m)
    {
        j = 0;
        while(j<n)
        {
            
            /*
             * Choose between specialized 4x4 code and general code
             */
            if( i+4<=m&&j+4<=n )
            {
                
                /*
                 * Specialized 4x4 code for [I..I+3]x[J..J+3] submatrix of C.
                 *
                 * This submatrix is calculated as sum of K rank-1 products,
                 * with operands cached in local variables in order to speed
                 * up operations with arrays.
                 */
                idxa0 = ia+i+0;
                idxa1 = ia+i+1;
                idxa2 = ia+i+2;
                idxa3 = ia+i+3;
                offsa = ja;
                idxb0 = ib+j+0;
                idxb1 = ib+j+1;
                idxb2 = ib+j+2;
                idxb3 = ib+j+3;
                offsb = jb;
                v00 = 0.0;
                v01 = 0.0;
                v02 = 0.0;
                v03 = 0.0;
                v10 = 0.0;
                v11 = 0.0;
                v12 = 0.0;
                v13 = 0.0;
                v20 = 0.0;
                v21 = 0.0;
                v22 = 0.0;
                v23 = 0.0;
                v30 = 0.0;
                v31 = 0.0;
                v32 = 0.0;
                v33 = 0.0;
                for(t=0; t<=k-1; t++)
                {
                    a0 = a->ptr.pp_double[idxa0][offsa];
                    a1 = a->ptr.pp_double[idxa1][offsa];
                    b0 = b->ptr.pp_double[idxb0][offsb];
                    b1 = b->ptr.pp_double[idxb1][offsb];
                    v00 = v00+a0*b0;
                    v01 = v01+a0*b1;
                    v10 = v10+a1*b0;
                    v11 = v11+a1*b1;
                    a2 = a->ptr.pp_double[idxa2][offsa];
                    a3 = a->ptr.pp_double[idxa3][offsa];
                    v20 = v20+a2*b0;
                    v21 = v21+a2*b1;
                    v30 = v30+a3*b0;
                    v31 = v31+a3*b1;
                    b2 = b->ptr.pp_double[idxb2][offsb];
                    b3 = b->ptr.pp_double[idxb3][offsb];
                    v22 = v22+a2*b2;
                    v23 = v23+a2*b3;
                    v32 = v32+a3*b2;
                    v33 = v33+a3*b3;
                    v02 = v02+a0*b2;
                    v03 = v03+a0*b3;
                    v12 = v12+a1*b2;
                    v13 = v13+a1*b3;
                    offsa = offsa+1;
                    offsb = offsb+1;
                }
                if( ae_fp_eq(beta,(double)(0)) )
                {
                    c->ptr.pp_double[ic+i+0][jc+j+0] = alpha*v00;
                    c->ptr.pp_double[ic+i+0][jc+j+1] = alpha*v01;
                    c->ptr.pp_double[ic+i+0][jc+j+2] = alpha*v02;
                    c->ptr.pp_double[ic+i+0][jc+j+3] = alpha*v03;
                    c->ptr.pp_double[ic+i+1][jc+j+0] = alpha*v10;
                    c->ptr.pp_double[ic+i+1][jc+j+1] = alpha*v11;
                    c->ptr.pp_double[ic+i+1][jc+j+2] = alpha*v12;
                    c->ptr.pp_double[ic+i+1][jc+j+3] = alpha*v13;
                    c->ptr.pp_double[ic+i+2][jc+j+0] = alpha*v20;
                    c->ptr.pp_double[ic+i+2][jc+j+1] = alpha*v21;
                    c->ptr.pp_double[ic+i+2][jc+j+2] = alpha*v22;
                    c->ptr.pp_double[ic+i+2][jc+j+3] = alpha*v23;
                    c->ptr.pp_double[ic+i+3][jc+j+0] = alpha*v30;
                    c->ptr.pp_double[ic+i+3][jc+j+1] = alpha*v31;
                    c->ptr.pp_double[ic+i+3][jc+j+2] = alpha*v32;
                    c->ptr.pp_double[ic+i+3][jc+j+3] = alpha*v33;
                }
                else
                {
                    c->ptr.pp_double[ic+i+0][jc+j+0] = beta*c->ptr.pp_double[ic+i+0][jc+j+0]+alpha*v00;
                    c->ptr.pp_double[ic+i+0][jc+j+1] = beta*c->ptr.pp_double[ic+i+0][jc+j+1]+alpha*v01;
                    c->ptr.pp_double[ic+i+0][jc+j+2] = beta*c->ptr.pp_double[ic+i+0][jc+j+2]+alpha*v02;
                    c->ptr.pp_double[ic+i+0][jc+j+3] = beta*c->ptr.pp_double[ic+i+0][jc+j+3]+alpha*v03;
                    c->ptr.pp_double[ic+i+1][jc+j+0] = beta*c->ptr.pp_double[ic+i+1][jc+j+0]+alpha*v10;
                    c->ptr.pp_double[ic+i+1][jc+j+1] = beta*c->ptr.pp_double[ic+i+1][jc+j+1]+alpha*v11;
                    c->ptr.pp_double[ic+i+1][jc+j+2] = beta*c->ptr.pp_double[ic+i+1][jc+j+2]+alpha*v12;
                    c->ptr.pp_double[ic+i+1][jc+j+3] = beta*c->ptr.pp_double[ic+i+1][jc+j+3]+alpha*v13;
                    c->ptr.pp_double[ic+i+2][jc+j+0] = beta*c->ptr.pp_double[ic+i+2][jc+j+0]+alpha*v20;
                    c->ptr.pp_double[ic+i+2][jc+j+1] = beta*c->ptr.pp_double[ic+i+2][jc+j+1]+alpha*v21;
                    c->ptr.pp_double[ic+i+2][jc+j+2] = beta*c->ptr.pp_double[ic+i+2][jc+j+2]+alpha*v22;
                    c->ptr.pp_double[ic+i+2][jc+j+3] = beta*c->ptr.pp_double[ic+i+2][jc+j+3]+alpha*v23;
                    c->ptr.pp_double[ic+i+3][jc+j+0] = beta*c->ptr.pp_double[ic+i+3][jc+j+0]+alpha*v30;
                    c->ptr.pp_double[ic+i+3][jc+j+1] = beta*c->ptr.pp_double[ic+i+3][jc+j+1]+alpha*v31;
                    c->ptr.pp_double[ic+i+3][jc+j+2] = beta*c->ptr.pp_double[ic+i+3][jc+j+2]+alpha*v32;
                    c->ptr.pp_double[ic+i+3][jc+j+3] = beta*c->ptr.pp_double[ic+i+3][jc+j+3]+alpha*v33;
                }
            }
            else
            {
                
                /*
                 * Determine submatrix [I0..I1]x[J0..J1] to process
                 */
                i0 = i;
                i1 = ae_minint(i+3, m-1, _state);
                j0 = j;
                j1 = ae_minint(j+3, n-1, _state);
                
                /*
                 * Process submatrix
                 */
                for(ik=i0; ik<=i1; ik++)
                {
                    for(jk=j0; jk<=j1; jk++)
                    {
                        if( k==0||ae_fp_eq(alpha,(double)(0)) )
                        {
                            v = (double)(0);
                        }
                        else
                        {
                            v = ae_v_dotproduct(&a->ptr.pp_double[ia+ik][ja], 1, &b->ptr.pp_double[ib+jk][jb], 1, ae_v_len(ja,ja+k-1));
                        }
                        if( ae_fp_eq(beta,(double)(0)) )
                        {
                            c->ptr.pp_double[ic+ik][jc+jk] = alpha*v;
                        }
                        else
                        {
                            c->ptr.pp_double[ic+ik][jc+jk] = beta*c->ptr.pp_double[ic+ik][jc+jk]+alpha*v;
                        }
                    }
                }
            }
            j = j+4;
        }
        i = i+4;
    }
}


/*************************************************************************
RMatrixGEMM kernel, basecase code for RMatrixGEMM, specialized for sitation
with OpTypeA=1 and OpTypeB=0.

Additional info:
* this function requires that Alpha<>0 (assertion is thrown otherwise)

INPUT PARAMETERS
    M       -   matrix size, M>0
    N       -   matrix size, N>0
    K       -   matrix size, K>0
    Alpha   -   coefficient
    A       -   matrix
    IA      -   submatrix offset
    JA      -   submatrix offset
    B       -   matrix
    IB      -   submatrix offset
    JB      -   submatrix offset
    Beta    -   coefficient
    C       -   PREALLOCATED output matrix
    IC      -   submatrix offset
    JC      -   submatrix offset

  -- ALGLIB routine --
     27.03.2013
     Bochkanov Sergey
*************************************************************************/
void rmatrixgemmk44v10(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     /* Real    */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double v;
    double v00;
    double v01;
    double v02;
    double v03;
    double v10;
    double v11;
    double v12;
    double v13;
    double v20;
    double v21;
    double v22;
    double v23;
    double v30;
    double v31;
    double v32;
    double v33;
    double a0;
    double a1;
    double a2;
    double a3;
    double b0;
    double b1;
    double b2;
    double b3;
    ae_int_t idxa0;
    ae_int_t idxa1;
    ae_int_t idxa2;
    ae_int_t idxa3;
    ae_int_t idxb0;
    ae_int_t idxb1;
    ae_int_t idxb2;
    ae_int_t idxb3;
    ae_int_t i0;
    ae_int_t i1;
    ae_int_t ik;
    ae_int_t j0;
    ae_int_t j1;
    ae_int_t jk;
    ae_int_t t;
    ae_int_t offsa;
    ae_int_t offsb;


    ae_assert(ae_fp_neq(alpha,(double)(0)), "RMatrixGEMMK44V00: internal error (Alpha=0)", _state);
    
    /*
     * if matrix size is zero
     */
    if( m==0||n==0 )
    {
        return;
    }
    
    /*
     * A'*B
     */
    i = 0;
    while(i<m)
    {
        j = 0;
        while(j<n)
        {
            
            /*
             * Choose between specialized 4x4 code and general code
             */
            if( i+4<=m&&j+4<=n )
            {
                
                /*
                 * Specialized 4x4 code for [I..I+3]x[J..J+3] submatrix of C.
                 *
                 * This submatrix is calculated as sum of K rank-1 products,
                 * with operands cached in local variables in order to speed
                 * up operations with arrays.
                 */
                idxa0 = ja+i+0;
                idxa1 = ja+i+1;
                idxa2 = ja+i+2;
                idxa3 = ja+i+3;
                offsa = ia;
                idxb0 = jb+j+0;
                idxb1 = jb+j+1;
                idxb2 = jb+j+2;
                idxb3 = jb+j+3;
                offsb = ib;
                v00 = 0.0;
                v01 = 0.0;
                v02 = 0.0;
                v03 = 0.0;
                v10 = 0.0;
                v11 = 0.0;
                v12 = 0.0;
                v13 = 0.0;
                v20 = 0.0;
                v21 = 0.0;
                v22 = 0.0;
                v23 = 0.0;
                v30 = 0.0;
                v31 = 0.0;
                v32 = 0.0;
                v33 = 0.0;
                for(t=0; t<=k-1; t++)
                {
                    a0 = a->ptr.pp_double[offsa][idxa0];
                    a1 = a->ptr.pp_double[offsa][idxa1];
                    b0 = b->ptr.pp_double[offsb][idxb0];
                    b1 = b->ptr.pp_double[offsb][idxb1];
                    v00 = v00+a0*b0;
                    v01 = v01+a0*b1;
                    v10 = v10+a1*b0;
                    v11 = v11+a1*b1;
                    a2 = a->ptr.pp_double[offsa][idxa2];
                    a3 = a->ptr.pp_double[offsa][idxa3];
                    v20 = v20+a2*b0;
                    v21 = v21+a2*b1;
                    v30 = v30+a3*b0;
                    v31 = v31+a3*b1;
                    b2 = b->ptr.pp_double[offsb][idxb2];
                    b3 = b->ptr.pp_double[offsb][idxb3];
                    v22 = v22+a2*b2;
                    v23 = v23+a2*b3;
                    v32 = v32+a3*b2;
                    v33 = v33+a3*b3;
                    v02 = v02+a0*b2;
                    v03 = v03+a0*b3;
                    v12 = v12+a1*b2;
                    v13 = v13+a1*b3;
                    offsa = offsa+1;
                    offsb = offsb+1;
                }
                if( ae_fp_eq(beta,(double)(0)) )
                {
                    c->ptr.pp_double[ic+i+0][jc+j+0] = alpha*v00;
                    c->ptr.pp_double[ic+i+0][jc+j+1] = alpha*v01;
                    c->ptr.pp_double[ic+i+0][jc+j+2] = alpha*v02;
                    c->ptr.pp_double[ic+i+0][jc+j+3] = alpha*v03;
                    c->ptr.pp_double[ic+i+1][jc+j+0] = alpha*v10;
                    c->ptr.pp_double[ic+i+1][jc+j+1] = alpha*v11;
                    c->ptr.pp_double[ic+i+1][jc+j+2] = alpha*v12;
                    c->ptr.pp_double[ic+i+1][jc+j+3] = alpha*v13;
                    c->ptr.pp_double[ic+i+2][jc+j+0] = alpha*v20;
                    c->ptr.pp_double[ic+i+2][jc+j+1] = alpha*v21;
                    c->ptr.pp_double[ic+i+2][jc+j+2] = alpha*v22;
                    c->ptr.pp_double[ic+i+2][jc+j+3] = alpha*v23;
                    c->ptr.pp_double[ic+i+3][jc+j+0] = alpha*v30;
                    c->ptr.pp_double[ic+i+3][jc+j+1] = alpha*v31;
                    c->ptr.pp_double[ic+i+3][jc+j+2] = alpha*v32;
                    c->ptr.pp_double[ic+i+3][jc+j+3] = alpha*v33;
                }
                else
                {
                    c->ptr.pp_double[ic+i+0][jc+j+0] = beta*c->ptr.pp_double[ic+i+0][jc+j+0]+alpha*v00;
                    c->ptr.pp_double[ic+i+0][jc+j+1] = beta*c->ptr.pp_double[ic+i+0][jc+j+1]+alpha*v01;
                    c->ptr.pp_double[ic+i+0][jc+j+2] = beta*c->ptr.pp_double[ic+i+0][jc+j+2]+alpha*v02;
                    c->ptr.pp_double[ic+i+0][jc+j+3] = beta*c->ptr.pp_double[ic+i+0][jc+j+3]+alpha*v03;
                    c->ptr.pp_double[ic+i+1][jc+j+0] = beta*c->ptr.pp_double[ic+i+1][jc+j+0]+alpha*v10;
                    c->ptr.pp_double[ic+i+1][jc+j+1] = beta*c->ptr.pp_double[ic+i+1][jc+j+1]+alpha*v11;
                    c->ptr.pp_double[ic+i+1][jc+j+2] = beta*c->ptr.pp_double[ic+i+1][jc+j+2]+alpha*v12;
                    c->ptr.pp_double[ic+i+1][jc+j+3] = beta*c->ptr.pp_double[ic+i+1][jc+j+3]+alpha*v13;
                    c->ptr.pp_double[ic+i+2][jc+j+0] = beta*c->ptr.pp_double[ic+i+2][jc+j+0]+alpha*v20;
                    c->ptr.pp_double[ic+i+2][jc+j+1] = beta*c->ptr.pp_double[ic+i+2][jc+j+1]+alpha*v21;
                    c->ptr.pp_double[ic+i+2][jc+j+2] = beta*c->ptr.pp_double[ic+i+2][jc+j+2]+alpha*v22;
                    c->ptr.pp_double[ic+i+2][jc+j+3] = beta*c->ptr.pp_double[ic+i+2][jc+j+3]+alpha*v23;
                    c->ptr.pp_double[ic+i+3][jc+j+0] = beta*c->ptr.pp_double[ic+i+3][jc+j+0]+alpha*v30;
                    c->ptr.pp_double[ic+i+3][jc+j+1] = beta*c->ptr.pp_double[ic+i+3][jc+j+1]+alpha*v31;
                    c->ptr.pp_double[ic+i+3][jc+j+2] = beta*c->ptr.pp_double[ic+i+3][jc+j+2]+alpha*v32;
                    c->ptr.pp_double[ic+i+3][jc+j+3] = beta*c->ptr.pp_double[ic+i+3][jc+j+3]+alpha*v33;
                }
            }
            else
            {
                
                /*
                 * Determine submatrix [I0..I1]x[J0..J1] to process
                 */
                i0 = i;
                i1 = ae_minint(i+3, m-1, _state);
                j0 = j;
                j1 = ae_minint(j+3, n-1, _state);
                
                /*
                 * Process submatrix
                 */
                for(ik=i0; ik<=i1; ik++)
                {
                    for(jk=j0; jk<=j1; jk++)
                    {
                        if( k==0||ae_fp_eq(alpha,(double)(0)) )
                        {
                            v = (double)(0);
                        }
                        else
                        {
                            v = 0.0;
                            v = ae_v_dotproduct(&a->ptr.pp_double[ia][ja+ik], a->stride, &b->ptr.pp_double[ib][jb+jk], b->stride, ae_v_len(ia,ia+k-1));
                        }
                        if( ae_fp_eq(beta,(double)(0)) )
                        {
                            c->ptr.pp_double[ic+ik][jc+jk] = alpha*v;
                        }
                        else
                        {
                            c->ptr.pp_double[ic+ik][jc+jk] = beta*c->ptr.pp_double[ic+ik][jc+jk]+alpha*v;
                        }
                    }
                }
            }
            j = j+4;
        }
        i = i+4;
    }
}


/*************************************************************************
RMatrixGEMM kernel, basecase code for RMatrixGEMM, specialized for sitation
with OpTypeA=1 and OpTypeB=1.

Additional info:
* this function requires that Alpha<>0 (assertion is thrown otherwise)

INPUT PARAMETERS
    M       -   matrix size, M>0
    N       -   matrix size, N>0
    K       -   matrix size, K>0
    Alpha   -   coefficient
    A       -   matrix
    IA      -   submatrix offset
    JA      -   submatrix offset
    B       -   matrix
    IB      -   submatrix offset
    JB      -   submatrix offset
    Beta    -   coefficient
    C       -   PREALLOCATED output matrix
    IC      -   submatrix offset
    JC      -   submatrix offset

  -- ALGLIB routine --
     27.03.2013
     Bochkanov Sergey
*************************************************************************/
void rmatrixgemmk44v11(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     /* Real    */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double v;
    double v00;
    double v01;
    double v02;
    double v03;
    double v10;
    double v11;
    double v12;
    double v13;
    double v20;
    double v21;
    double v22;
    double v23;
    double v30;
    double v31;
    double v32;
    double v33;
    double a0;
    double a1;
    double a2;
    double a3;
    double b0;
    double b1;
    double b2;
    double b3;
    ae_int_t idxa0;
    ae_int_t idxa1;
    ae_int_t idxa2;
    ae_int_t idxa3;
    ae_int_t idxb0;
    ae_int_t idxb1;
    ae_int_t idxb2;
    ae_int_t idxb3;
    ae_int_t i0;
    ae_int_t i1;
    ae_int_t ik;
    ae_int_t j0;
    ae_int_t j1;
    ae_int_t jk;
    ae_int_t t;
    ae_int_t offsa;
    ae_int_t offsb;


    ae_assert(ae_fp_neq(alpha,(double)(0)), "RMatrixGEMMK44V00: internal error (Alpha=0)", _state);
    
    /*
     * if matrix size is zero
     */
    if( m==0||n==0 )
    {
        return;
    }
    
    /*
     * A'*B'
     */
    i = 0;
    while(i<m)
    {
        j = 0;
        while(j<n)
        {
            
            /*
             * Choose between specialized 4x4 code and general code
             */
            if( i+4<=m&&j+4<=n )
            {
                
                /*
                 * Specialized 4x4 code for [I..I+3]x[J..J+3] submatrix of C.
                 *
                 * This submatrix is calculated as sum of K rank-1 products,
                 * with operands cached in local variables in order to speed
                 * up operations with arrays.
                 */
                idxa0 = ja+i+0;
                idxa1 = ja+i+1;
                idxa2 = ja+i+2;
                idxa3 = ja+i+3;
                offsa = ia;
                idxb0 = ib+j+0;
                idxb1 = ib+j+1;
                idxb2 = ib+j+2;
                idxb3 = ib+j+3;
                offsb = jb;
                v00 = 0.0;
                v01 = 0.0;
                v02 = 0.0;
                v03 = 0.0;
                v10 = 0.0;
                v11 = 0.0;
                v12 = 0.0;
                v13 = 0.0;
                v20 = 0.0;
                v21 = 0.0;
                v22 = 0.0;
                v23 = 0.0;
                v30 = 0.0;
                v31 = 0.0;
                v32 = 0.0;
                v33 = 0.0;
                for(t=0; t<=k-1; t++)
                {
                    a0 = a->ptr.pp_double[offsa][idxa0];
                    a1 = a->ptr.pp_double[offsa][idxa1];
                    b0 = b->ptr.pp_double[idxb0][offsb];
                    b1 = b->ptr.pp_double[idxb1][offsb];
                    v00 = v00+a0*b0;
                    v01 = v01+a0*b1;
                    v10 = v10+a1*b0;
                    v11 = v11+a1*b1;
                    a2 = a->ptr.pp_double[offsa][idxa2];
                    a3 = a->ptr.pp_double[offsa][idxa3];
                    v20 = v20+a2*b0;
                    v21 = v21+a2*b1;
                    v30 = v30+a3*b0;
                    v31 = v31+a3*b1;
                    b2 = b->ptr.pp_double[idxb2][offsb];
                    b3 = b->ptr.pp_double[idxb3][offsb];
                    v22 = v22+a2*b2;
                    v23 = v23+a2*b3;
                    v32 = v32+a3*b2;
                    v33 = v33+a3*b3;
                    v02 = v02+a0*b2;
                    v03 = v03+a0*b3;
                    v12 = v12+a1*b2;
                    v13 = v13+a1*b3;
                    offsa = offsa+1;
                    offsb = offsb+1;
                }
                if( ae_fp_eq(beta,(double)(0)) )
                {
                    c->ptr.pp_double[ic+i+0][jc+j+0] = alpha*v00;
                    c->ptr.pp_double[ic+i+0][jc+j+1] = alpha*v01;
                    c->ptr.pp_double[ic+i+0][jc+j+2] = alpha*v02;
                    c->ptr.pp_double[ic+i+0][jc+j+3] = alpha*v03;
                    c->ptr.pp_double[ic+i+1][jc+j+0] = alpha*v10;
                    c->ptr.pp_double[ic+i+1][jc+j+1] = alpha*v11;
                    c->ptr.pp_double[ic+i+1][jc+j+2] = alpha*v12;
                    c->ptr.pp_double[ic+i+1][jc+j+3] = alpha*v13;
                    c->ptr.pp_double[ic+i+2][jc+j+0] = alpha*v20;
                    c->ptr.pp_double[ic+i+2][jc+j+1] = alpha*v21;
                    c->ptr.pp_double[ic+i+2][jc+j+2] = alpha*v22;
                    c->ptr.pp_double[ic+i+2][jc+j+3] = alpha*v23;
                    c->ptr.pp_double[ic+i+3][jc+j+0] = alpha*v30;
                    c->ptr.pp_double[ic+i+3][jc+j+1] = alpha*v31;
                    c->ptr.pp_double[ic+i+3][jc+j+2] = alpha*v32;
                    c->ptr.pp_double[ic+i+3][jc+j+3] = alpha*v33;
                }
                else
                {
                    c->ptr.pp_double[ic+i+0][jc+j+0] = beta*c->ptr.pp_double[ic+i+0][jc+j+0]+alpha*v00;
                    c->ptr.pp_double[ic+i+0][jc+j+1] = beta*c->ptr.pp_double[ic+i+0][jc+j+1]+alpha*v01;
                    c->ptr.pp_double[ic+i+0][jc+j+2] = beta*c->ptr.pp_double[ic+i+0][jc+j+2]+alpha*v02;
                    c->ptr.pp_double[ic+i+0][jc+j+3] = beta*c->ptr.pp_double[ic+i+0][jc+j+3]+alpha*v03;
                    c->ptr.pp_double[ic+i+1][jc+j+0] = beta*c->ptr.pp_double[ic+i+1][jc+j+0]+alpha*v10;
                    c->ptr.pp_double[ic+i+1][jc+j+1] = beta*c->ptr.pp_double[ic+i+1][jc+j+1]+alpha*v11;
                    c->ptr.pp_double[ic+i+1][jc+j+2] = beta*c->ptr.pp_double[ic+i+1][jc+j+2]+alpha*v12;
                    c->ptr.pp_double[ic+i+1][jc+j+3] = beta*c->ptr.pp_double[ic+i+1][jc+j+3]+alpha*v13;
                    c->ptr.pp_double[ic+i+2][jc+j+0] = beta*c->ptr.pp_double[ic+i+2][jc+j+0]+alpha*v20;
                    c->ptr.pp_double[ic+i+2][jc+j+1] = beta*c->ptr.pp_double[ic+i+2][jc+j+1]+alpha*v21;
                    c->ptr.pp_double[ic+i+2][jc+j+2] = beta*c->ptr.pp_double[ic+i+2][jc+j+2]+alpha*v22;
                    c->ptr.pp_double[ic+i+2][jc+j+3] = beta*c->ptr.pp_double[ic+i+2][jc+j+3]+alpha*v23;
                    c->ptr.pp_double[ic+i+3][jc+j+0] = beta*c->ptr.pp_double[ic+i+3][jc+j+0]+alpha*v30;
                    c->ptr.pp_double[ic+i+3][jc+j+1] = beta*c->ptr.pp_double[ic+i+3][jc+j+1]+alpha*v31;
                    c->ptr.pp_double[ic+i+3][jc+j+2] = beta*c->ptr.pp_double[ic+i+3][jc+j+2]+alpha*v32;
                    c->ptr.pp_double[ic+i+3][jc+j+3] = beta*c->ptr.pp_double[ic+i+3][jc+j+3]+alpha*v33;
                }
            }
            else
            {
                
                /*
                 * Determine submatrix [I0..I1]x[J0..J1] to process
                 */
                i0 = i;
                i1 = ae_minint(i+3, m-1, _state);
                j0 = j;
                j1 = ae_minint(j+3, n-1, _state);
                
                /*
                 * Process submatrix
                 */
                for(ik=i0; ik<=i1; ik++)
                {
                    for(jk=j0; jk<=j1; jk++)
                    {
                        if( k==0||ae_fp_eq(alpha,(double)(0)) )
                        {
                            v = (double)(0);
                        }
                        else
                        {
                            v = 0.0;
                            v = ae_v_dotproduct(&a->ptr.pp_double[ia][ja+ik], a->stride, &b->ptr.pp_double[ib+jk][jb], 1, ae_v_len(ia,ia+k-1));
                        }
                        if( ae_fp_eq(beta,(double)(0)) )
                        {
                            c->ptr.pp_double[ic+ik][jc+jk] = alpha*v;
                        }
                        else
                        {
                            c->ptr.pp_double[ic+ik][jc+jk] = beta*c->ptr.pp_double[ic+ik][jc+jk]+alpha*v;
                        }
                    }
                }
            }
            j = j+4;
        }
        i = i+4;
    }
}




/*************************************************************************
MKL-based kernel

  -- ALGLIB routine --
     01.10.2013
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixsyrkmkl(ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_bool isupper,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixsyrkmkl(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}


/*************************************************************************
MKL-based kernel

  -- ALGLIB routine --
     01.10.2013
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixherkmkl(ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     double beta,
     /* Complex */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_bool isupper,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixherkmkl(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}


/*************************************************************************
MKL-based kernel

  -- ALGLIB routine --
     01.10.2013
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixgemmmkl(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     /* Real    */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     ae_int_t optypeb,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixgemmmkl(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}


/*************************************************************************
MKL-based kernel

  -- ALGLIB routine --
     16.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixgemmmkl(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     ae_complex alpha,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     /* Complex */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     ae_int_t optypeb,
     ae_complex beta,
     /* Complex */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixgemmmkl(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}


/*************************************************************************
MKL-based kernel

  -- ALGLIB routine --
     16.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixlefttrsmmkl(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Complex */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixlefttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
MKL-based kernel

  -- ALGLIB routine --
     16.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixrighttrsmmkl(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Complex */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixrighttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
MKL-based kernel

  -- ALGLIB routine --
     16.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixlefttrsmmkl(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Real    */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixlefttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
MKL-based kernel

  -- ALGLIB routine --
     16.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixrighttrsmmkl(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Real    */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixrighttrsmmkl(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
MKL-based kernel.

NOTE:

if function returned False, CholResult is NOT modified. Not ever referenced!
if function returned True, CholResult is set to status of Cholesky decomposition
(True on succeess).

  -- ALGLIB routine --
     16.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool spdmatrixcholeskymkl(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t n,
     ae_bool isupper,
     ae_bool* cholresult,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_spdmatrixcholeskymkl(a, offs, n, isupper, cholresult);
#endif
}


/*************************************************************************
MKL-based kernel.

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixplumkl(/* Real    */ ae_matrix* a,
     ae_int_t offs,
     ae_int_t m,
     ae_int_t n,
     /* Integer */ ae_vector* pivots,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixplumkl(a, offs, m, n, pivots);
#endif
}


/*************************************************************************
MKL-based kernel.

NOTE: this function needs preallocated output/temporary arrays.
      D and E must be at least max(M,N)-wide.

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixbdmkl(/* Real    */ ae_matrix* a,
     ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     /* Real    */ ae_vector* tauq,
     /* Real    */ ae_vector* taup,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixbdmkl(a, m, n, d, e, tauq, taup);
#endif
}


/*************************************************************************
MKL-based kernel.

If ByQ is True,  TauP is not used (can be empty array).
If ByQ is False, TauQ is not used (can be empty array).

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixbdmultiplybymkl(/* Real    */ ae_matrix* qp,
     ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_vector* tauq,
     /* Real    */ ae_vector* taup,
     /* Real    */ ae_matrix* z,
     ae_int_t zrows,
     ae_int_t zcolumns,
     ae_bool byq,
     ae_bool fromtheright,
     ae_bool dotranspose,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixbdmultiplybymkl(qp, m, n, tauq, taup, z, zrows, zcolumns, byq, fromtheright, dotranspose);
#endif
}


/*************************************************************************
MKL-based kernel.

NOTE: Tau must be preallocated array with at least N-1 elements.

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixhessenbergmkl(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* tau,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixhessenbergmkl(a, n, tau);
#endif
}


/*************************************************************************
MKL-based kernel.

NOTE: Q must be preallocated N*N array

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixhessenbergunpackqmkl(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* tau,
     /* Real    */ ae_matrix* q,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixhessenbergunpackqmkl(a, n, tau, q);
#endif
}


/*************************************************************************
MKL-based kernel.

NOTE: Tau, D, E must be preallocated arrays;
      length(E)=length(Tau)=N-1 (or larger)
      length(D)=N (or larger)

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool smatrixtdmkl(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* tau,
     /* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_smatrixtdmkl(a, n, isupper, tau, d, e);
#endif
}


/*************************************************************************
MKL-based kernel.

NOTE: Q must be preallocated N*N array

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool smatrixtdunpackqmkl(/* Real    */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_vector* tau,
     /* Real    */ ae_matrix* q,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_smatrixtdunpackqmkl(a, n, isupper, tau, q);
#endif
}


/*************************************************************************
MKL-based kernel.

NOTE: Tau, D, E must be preallocated arrays;
      length(E)=length(Tau)=N-1 (or larger)
      length(D)=N (or larger)

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool hmatrixtdmkl(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* tau,
     /* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_hmatrixtdmkl(a, n, isupper, tau, d, e);
#endif
}


/*************************************************************************
MKL-based kernel.

NOTE: Q must be preallocated N*N array

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool hmatrixtdunpackqmkl(/* Complex */ ae_matrix* a,
     ae_int_t n,
     ae_bool isupper,
     /* Complex */ ae_vector* tau,
     /* Complex */ ae_matrix* q,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_hmatrixtdunpackqmkl(a, n, isupper, tau, q);
#endif
}


/*************************************************************************
MKL-based kernel.

Returns True if MKL was present and handled request (MKL  completion  code
is returned as separate output parameter).

D and E are pre-allocated arrays with length N (both of them!). On output,
D constraints singular values, and E is destroyed.

SVDResult is modified if and only if MKL is present.

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixbdsvdmkl(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     ae_bool isupper,
     /* Real    */ ae_matrix* u,
     ae_int_t nru,
     /* Real    */ ae_matrix* c,
     ae_int_t ncc,
     /* Real    */ ae_matrix* vt,
     ae_int_t ncvt,
     ae_bool* svdresult,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixbdsvdmkl(d, e, n, isupper, u, nru, c, ncc, vt, ncvt, svdresult);
#endif
}


/*************************************************************************
MKL-based DHSEQR kernel.

Returns True if MKL was present and handled request.

WR and WI are pre-allocated arrays with length N.
Z is pre-allocated array[N,N].

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixinternalschurdecompositionmkl(/* Real    */ ae_matrix* h,
     ae_int_t n,
     ae_int_t tneeded,
     ae_int_t zneeded,
     /* Real    */ ae_vector* wr,
     /* Real    */ ae_vector* wi,
     /* Real    */ ae_matrix* z,
     ae_int_t* info,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixinternalschurdecompositionmkl(h, n, tneeded, zneeded, wr, wi, z, info);
#endif
}


/*************************************************************************
MKL-based DTREVC kernel.

Returns True if MKL was present and handled request.

NOTE: this function does NOT support HOWMNY=3!!!!

VL and VR are pre-allocated arrays with length N*N, if required. If particalar
variables is not required, it can be dummy (empty) array.

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixinternaltrevcmkl(/* Real    */ ae_matrix* t,
     ae_int_t n,
     ae_int_t side,
     ae_int_t howmny,
     /* Real    */ ae_matrix* vl,
     /* Real    */ ae_matrix* vr,
     ae_int_t* m,
     ae_int_t* info,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixinternaltrevcmkl(t, n, side, howmny, vl, vr, m, info);
#endif
}


/*************************************************************************
MKL-based kernel.

Returns True if MKL was present and handled request (MKL  completion  code
is returned as separate output parameter).

D and E are pre-allocated arrays with length N (both of them!). On output,
D constraints eigenvalues, and E is destroyed.

Z is preallocated array[N,N] for ZNeeded<>0; ignored for ZNeeded=0.

EVDResult is modified if and only if MKL is present.

  -- ALGLIB routine --
     20.10.2014
     Bochkanov Sergey
*************************************************************************/
ae_bool smatrixtdevdmkl(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* e,
     ae_int_t n,
     ae_int_t zneeded,
     /* Real    */ ae_matrix* z,
     ae_bool* evdresult,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_MKL
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_smatrixtdevdmkl(d, e, n, zneeded, z, evdresult);
#endif
}




double vectornorm2(/* Real    */ ae_vector* x,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t ix;
    double absxi;
    double scl;
    double ssq;
    double result;


    n = i2-i1+1;
    if( n<1 )
    {
        result = (double)(0);
        return result;
    }
    if( n==1 )
    {
        result = ae_fabs(x->ptr.p_double[i1], _state);
        return result;
    }
    scl = (double)(0);
    ssq = (double)(1);
    for(ix=i1; ix<=i2; ix++)
    {
        if( ae_fp_neq(x->ptr.p_double[ix],(double)(0)) )
        {
            absxi = ae_fabs(x->ptr.p_double[ix], _state);
            if( ae_fp_less(scl,absxi) )
            {
                ssq = 1+ssq*ae_sqr(scl/absxi, _state);
                scl = absxi;
            }
            else
            {
                ssq = ssq+ae_sqr(absxi/scl, _state);
            }
        }
    }
    result = scl*ae_sqrt(ssq, _state);
    return result;
}


ae_int_t vectoridxabsmax(/* Real    */ ae_vector* x,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t result;


    result = i1;
    for(i=i1+1; i<=i2; i++)
    {
        if( ae_fp_greater(ae_fabs(x->ptr.p_double[i], _state),ae_fabs(x->ptr.p_double[result], _state)) )
        {
            result = i;
        }
    }
    return result;
}


ae_int_t columnidxabsmax(/* Real    */ ae_matrix* x,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t result;


    result = i1;
    for(i=i1+1; i<=i2; i++)
    {
        if( ae_fp_greater(ae_fabs(x->ptr.pp_double[i][j], _state),ae_fabs(x->ptr.pp_double[result][j], _state)) )
        {
            result = i;
        }
    }
    return result;
}


ae_int_t rowidxabsmax(/* Real    */ ae_matrix* x,
     ae_int_t j1,
     ae_int_t j2,
     ae_int_t i,
     ae_state *_state)
{
    ae_int_t j;
    ae_int_t result;


    result = j1;
    for(j=j1+1; j<=j2; j++)
    {
        if( ae_fp_greater(ae_fabs(x->ptr.pp_double[i][j], _state),ae_fabs(x->ptr.pp_double[i][result], _state)) )
        {
            result = j;
        }
    }
    return result;
}


double upperhessenberg1norm(/* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j1,
     ae_int_t j2,
     /* Real    */ ae_vector* work,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double result;


    ae_assert(i2-i1==j2-j1, "UpperHessenberg1Norm: I2-I1<>J2-J1!", _state);
    for(j=j1; j<=j2; j++)
    {
        work->ptr.p_double[j] = (double)(0);
    }
    for(i=i1; i<=i2; i++)
    {
        for(j=ae_maxint(j1, j1+i-i1-1, _state); j<=j2; j++)
        {
            work->ptr.p_double[j] = work->ptr.p_double[j]+ae_fabs(a->ptr.pp_double[i][j], _state);
        }
    }
    result = (double)(0);
    for(j=j1; j<=j2; j++)
    {
        result = ae_maxreal(result, work->ptr.p_double[j], _state);
    }
    return result;
}


void copymatrix(/* Real    */ ae_matrix* a,
     ae_int_t is1,
     ae_int_t is2,
     ae_int_t js1,
     ae_int_t js2,
     /* Real    */ ae_matrix* b,
     ae_int_t id1,
     ae_int_t id2,
     ae_int_t jd1,
     ae_int_t jd2,
     ae_state *_state)
{
    ae_int_t isrc;
    ae_int_t idst;


    if( is1>is2||js1>js2 )
    {
        return;
    }
    ae_assert(is2-is1==id2-id1, "CopyMatrix: different sizes!", _state);
    ae_assert(js2-js1==jd2-jd1, "CopyMatrix: different sizes!", _state);
    for(isrc=is1; isrc<=is2; isrc++)
    {
        idst = isrc-is1+id1;
        ae_v_move(&b->ptr.pp_double[idst][jd1], 1, &a->ptr.pp_double[isrc][js1], 1, ae_v_len(jd1,jd2));
    }
}


void inplacetranspose(/* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j1,
     ae_int_t j2,
     /* Real    */ ae_vector* work,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t ips;
    ae_int_t jps;
    ae_int_t l;


    if( i1>i2||j1>j2 )
    {
        return;
    }
    ae_assert(i1-i2==j1-j2, "InplaceTranspose error: incorrect array size!", _state);
    for(i=i1; i<=i2-1; i++)
    {
        j = j1+i-i1;
        ips = i+1;
        jps = j1+ips-i1;
        l = i2-i;
        ae_v_move(&work->ptr.p_double[1], 1, &a->ptr.pp_double[ips][j], a->stride, ae_v_len(1,l));
        ae_v_move(&a->ptr.pp_double[ips][j], a->stride, &a->ptr.pp_double[i][jps], 1, ae_v_len(ips,i2));
        ae_v_move(&a->ptr.pp_double[i][jps], 1, &work->ptr.p_double[1], 1, ae_v_len(jps,j2));
    }
}


void copyandtranspose(/* Real    */ ae_matrix* a,
     ae_int_t is1,
     ae_int_t is2,
     ae_int_t js1,
     ae_int_t js2,
     /* Real    */ ae_matrix* b,
     ae_int_t id1,
     ae_int_t id2,
     ae_int_t jd1,
     ae_int_t jd2,
     ae_state *_state)
{
    ae_int_t isrc;
    ae_int_t jdst;


    if( is1>is2||js1>js2 )
    {
        return;
    }
    ae_assert(is2-is1==jd2-jd1, "CopyAndTranspose: different sizes!", _state);
    ae_assert(js2-js1==id2-id1, "CopyAndTranspose: different sizes!", _state);
    for(isrc=is1; isrc<=is2; isrc++)
    {
        jdst = isrc-is1+jd1;
        ae_v_move(&b->ptr.pp_double[id1][jdst], b->stride, &a->ptr.pp_double[isrc][js1], 1, ae_v_len(id1,id2));
    }
}


void matrixvectormultiply(/* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j1,
     ae_int_t j2,
     ae_bool trans,
     /* Real    */ ae_vector* x,
     ae_int_t ix1,
     ae_int_t ix2,
     double alpha,
     /* Real    */ ae_vector* y,
     ae_int_t iy1,
     ae_int_t iy2,
     double beta,
     ae_state *_state)
{
    ae_int_t i;
    double v;


    if( !trans )
    {
        
        /*
         * y := alpha*A*x + beta*y;
         */
        if( i1>i2||j1>j2 )
        {
            return;
        }
        ae_assert(j2-j1==ix2-ix1, "MatrixVectorMultiply: A and X dont match!", _state);
        ae_assert(i2-i1==iy2-iy1, "MatrixVectorMultiply: A and Y dont match!", _state);
        
        /*
         * beta*y
         */
        if( ae_fp_eq(beta,(double)(0)) )
        {
            for(i=iy1; i<=iy2; i++)
            {
                y->ptr.p_double[i] = (double)(0);
            }
        }
        else
        {
            ae_v_muld(&y->ptr.p_double[iy1], 1, ae_v_len(iy1,iy2), beta);
        }
        
        /*
         * alpha*A*x
         */
        for(i=i1; i<=i2; i++)
        {
            v = ae_v_dotproduct(&a->ptr.pp_double[i][j1], 1, &x->ptr.p_double[ix1], 1, ae_v_len(j1,j2));
            y->ptr.p_double[iy1+i-i1] = y->ptr.p_double[iy1+i-i1]+alpha*v;
        }
    }
    else
    {
        
        /*
         * y := alpha*A'*x + beta*y;
         */
        if( i1>i2||j1>j2 )
        {
            return;
        }
        ae_assert(i2-i1==ix2-ix1, "MatrixVectorMultiply: A and X dont match!", _state);
        ae_assert(j2-j1==iy2-iy1, "MatrixVectorMultiply: A and Y dont match!", _state);
        
        /*
         * beta*y
         */
        if( ae_fp_eq(beta,(double)(0)) )
        {
            for(i=iy1; i<=iy2; i++)
            {
                y->ptr.p_double[i] = (double)(0);
            }
        }
        else
        {
            ae_v_muld(&y->ptr.p_double[iy1], 1, ae_v_len(iy1,iy2), beta);
        }
        
        /*
         * alpha*A'*x
         */
        for(i=i1; i<=i2; i++)
        {
            v = alpha*x->ptr.p_double[ix1+i-i1];
            ae_v_addd(&y->ptr.p_double[iy1], 1, &a->ptr.pp_double[i][j1], 1, ae_v_len(iy1,iy2), v);
        }
    }
}


double pythag2(double x, double y, ae_state *_state)
{
    double w;
    double xabs;
    double yabs;
    double z;
    double result;


    xabs = ae_fabs(x, _state);
    yabs = ae_fabs(y, _state);
    w = ae_maxreal(xabs, yabs, _state);
    z = ae_minreal(xabs, yabs, _state);
    if( ae_fp_eq(z,(double)(0)) )
    {
        result = w;
    }
    else
    {
        result = w*ae_sqrt(1+ae_sqr(z/w, _state), _state);
    }
    return result;
}


void matrixmatrixmultiply(/* Real    */ ae_matrix* a,
     ae_int_t ai1,
     ae_int_t ai2,
     ae_int_t aj1,
     ae_int_t aj2,
     ae_bool transa,
     /* Real    */ ae_matrix* b,
     ae_int_t bi1,
     ae_int_t bi2,
     ae_int_t bj1,
     ae_int_t bj2,
     ae_bool transb,
     double alpha,
     /* Real    */ ae_matrix* c,
     ae_int_t ci1,
     ae_int_t ci2,
     ae_int_t cj1,
     ae_int_t cj2,
     double beta,
     /* Real    */ ae_vector* work,
     ae_state *_state)
{
    ae_int_t arows;
    ae_int_t acols;
    ae_int_t brows;
    ae_int_t bcols;
    ae_int_t crows;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t l;
    ae_int_t r;
    double v;


    
    /*
     * Setup
     */
    if( !transa )
    {
        arows = ai2-ai1+1;
        acols = aj2-aj1+1;
    }
    else
    {
        arows = aj2-aj1+1;
        acols = ai2-ai1+1;
    }
    if( !transb )
    {
        brows = bi2-bi1+1;
        bcols = bj2-bj1+1;
    }
    else
    {
        brows = bj2-bj1+1;
        bcols = bi2-bi1+1;
    }
    ae_assert(acols==brows, "MatrixMatrixMultiply: incorrect matrix sizes!", _state);
    if( ((arows<=0||acols<=0)||brows<=0)||bcols<=0 )
    {
        return;
    }
    crows = arows;
    
    /*
     * Test WORK
     */
    i = ae_maxint(arows, acols, _state);
    i = ae_maxint(brows, i, _state);
    i = ae_maxint(i, bcols, _state);
    work->ptr.p_double[1] = (double)(0);
    work->ptr.p_double[i] = (double)(0);
    
    /*
     * Prepare C
     */
    if( ae_fp_eq(beta,(double)(0)) )
    {
        for(i=ci1; i<=ci2; i++)
        {
            for(j=cj1; j<=cj2; j++)
            {
                c->ptr.pp_double[i][j] = (double)(0);
            }
        }
    }
    else
    {
        for(i=ci1; i<=ci2; i++)
        {
            ae_v_muld(&c->ptr.pp_double[i][cj1], 1, ae_v_len(cj1,cj2), beta);
        }
    }
    
    /*
     * A*B
     */
    if( !transa&&!transb )
    {
        for(l=ai1; l<=ai2; l++)
        {
            for(r=bi1; r<=bi2; r++)
            {
                v = alpha*a->ptr.pp_double[l][aj1+r-bi1];
                k = ci1+l-ai1;
                ae_v_addd(&c->ptr.pp_double[k][cj1], 1, &b->ptr.pp_double[r][bj1], 1, ae_v_len(cj1,cj2), v);
            }
        }
        return;
    }
    
    /*
     * A*B'
     */
    if( !transa&&transb )
    {
        if( arows*acols<brows*bcols )
        {
            for(r=bi1; r<=bi2; r++)
            {
                for(l=ai1; l<=ai2; l++)
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[l][aj1], 1, &b->ptr.pp_double[r][bj1], 1, ae_v_len(aj1,aj2));
                    c->ptr.pp_double[ci1+l-ai1][cj1+r-bi1] = c->ptr.pp_double[ci1+l-ai1][cj1+r-bi1]+alpha*v;
                }
            }
            return;
        }
        else
        {
            for(l=ai1; l<=ai2; l++)
            {
                for(r=bi1; r<=bi2; r++)
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[l][aj1], 1, &b->ptr.pp_double[r][bj1], 1, ae_v_len(aj1,aj2));
                    c->ptr.pp_double[ci1+l-ai1][cj1+r-bi1] = c->ptr.pp_double[ci1+l-ai1][cj1+r-bi1]+alpha*v;
                }
            }
            return;
        }
    }
    
    /*
     * A'*B
     */
    if( transa&&!transb )
    {
        for(l=aj1; l<=aj2; l++)
        {
            for(r=bi1; r<=bi2; r++)
            {
                v = alpha*a->ptr.pp_double[ai1+r-bi1][l];
                k = ci1+l-aj1;
                ae_v_addd(&c->ptr.pp_double[k][cj1], 1, &b->ptr.pp_double[r][bj1], 1, ae_v_len(cj1,cj2), v);
            }
        }
        return;
    }
    
    /*
     * A'*B'
     */
    if( transa&&transb )
    {
        if( arows*acols<brows*bcols )
        {
            for(r=bi1; r<=bi2; r++)
            {
                k = cj1+r-bi1;
                for(i=1; i<=crows; i++)
                {
                    work->ptr.p_double[i] = 0.0;
                }
                for(l=ai1; l<=ai2; l++)
                {
                    v = alpha*b->ptr.pp_double[r][bj1+l-ai1];
                    ae_v_addd(&work->ptr.p_double[1], 1, &a->ptr.pp_double[l][aj1], 1, ae_v_len(1,crows), v);
                }
                ae_v_add(&c->ptr.pp_double[ci1][k], c->stride, &work->ptr.p_double[1], 1, ae_v_len(ci1,ci2));
            }
            return;
        }
        else
        {
            for(l=aj1; l<=aj2; l++)
            {
                k = ai2-ai1+1;
                ae_v_move(&work->ptr.p_double[1], 1, &a->ptr.pp_double[ai1][l], a->stride, ae_v_len(1,k));
                for(r=bi1; r<=bi2; r++)
                {
                    v = ae_v_dotproduct(&work->ptr.p_double[1], 1, &b->ptr.pp_double[r][bj1], 1, ae_v_len(1,k));
                    c->ptr.pp_double[ci1+l-aj1][cj1+r-bi1] = c->ptr.pp_double[ci1+l-aj1][cj1+r-bi1]+alpha*v;
                }
            }
            return;
        }
    }
}




void hermitianmatrixvectormultiply(/* Complex */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Complex */ ae_vector* x,
     ae_complex alpha,
     /* Complex */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ba1;
    ae_int_t by1;
    ae_int_t by2;
    ae_int_t bx1;
    ae_int_t bx2;
    ae_int_t n;
    ae_complex v;


    n = i2-i1+1;
    if( n<=0 )
    {
        return;
    }
    
    /*
     * Let A = L + D + U, where
     *  L is strictly lower triangular (main diagonal is zero)
     *  D is diagonal
     *  U is strictly upper triangular (main diagonal is zero)
     *
     * A*x = L*x + D*x + U*x
     *
     * Calculate D*x first
     */
    for(i=i1; i<=i2; i++)
    {
        y->ptr.p_complex[i-i1+1] = ae_c_mul(a->ptr.pp_complex[i][i],x->ptr.p_complex[i-i1+1]);
    }
    
    /*
     * Add L*x + U*x
     */
    if( isupper )
    {
        for(i=i1; i<=i2-1; i++)
        {
            
            /*
             * Add L*x to the result
             */
            v = x->ptr.p_complex[i-i1+1];
            by1 = i-i1+2;
            by2 = n;
            ba1 = i+1;
            ae_v_caddc(&y->ptr.p_complex[by1], 1, &a->ptr.pp_complex[i][ba1], 1, "Conj", ae_v_len(by1,by2), v);
            
            /*
             * Add U*x to the result
             */
            bx1 = i-i1+2;
            bx2 = n;
            ba1 = i+1;
            v = ae_v_cdotproduct(&x->ptr.p_complex[bx1], 1, "N", &a->ptr.pp_complex[i][ba1], 1, "N", ae_v_len(bx1,bx2));
            y->ptr.p_complex[i-i1+1] = ae_c_add(y->ptr.p_complex[i-i1+1],v);
        }
    }
    else
    {
        for(i=i1+1; i<=i2; i++)
        {
            
            /*
             * Add L*x to the result
             */
            bx1 = 1;
            bx2 = i-i1;
            ba1 = i1;
            v = ae_v_cdotproduct(&x->ptr.p_complex[bx1], 1, "N", &a->ptr.pp_complex[i][ba1], 1, "N", ae_v_len(bx1,bx2));
            y->ptr.p_complex[i-i1+1] = ae_c_add(y->ptr.p_complex[i-i1+1],v);
            
            /*
             * Add U*x to the result
             */
            v = x->ptr.p_complex[i-i1+1];
            by1 = 1;
            by2 = i-i1;
            ba1 = i1;
            ae_v_caddc(&y->ptr.p_complex[by1], 1, &a->ptr.pp_complex[i][ba1], 1, "Conj", ae_v_len(by1,by2), v);
        }
    }
    ae_v_cmulc(&y->ptr.p_complex[1], 1, ae_v_len(1,n), alpha);
}


void hermitianrank2update(/* Complex */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Complex */ ae_vector* x,
     /* Complex */ ae_vector* y,
     /* Complex */ ae_vector* t,
     ae_complex alpha,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t tp1;
    ae_int_t tp2;
    ae_complex v;


    if( isupper )
    {
        for(i=i1; i<=i2; i++)
        {
            tp1 = i+1-i1;
            tp2 = i2-i1+1;
            v = ae_c_mul(alpha,x->ptr.p_complex[i+1-i1]);
            ae_v_cmovec(&t->ptr.p_complex[tp1], 1, &y->ptr.p_complex[tp1], 1, "Conj", ae_v_len(tp1,tp2), v);
            v = ae_c_mul(ae_c_conj(alpha, _state),y->ptr.p_complex[i+1-i1]);
            ae_v_caddc(&t->ptr.p_complex[tp1], 1, &x->ptr.p_complex[tp1], 1, "Conj", ae_v_len(tp1,tp2), v);
            ae_v_cadd(&a->ptr.pp_complex[i][i], 1, &t->ptr.p_complex[tp1], 1, "N", ae_v_len(i,i2));
        }
    }
    else
    {
        for(i=i1; i<=i2; i++)
        {
            tp1 = 1;
            tp2 = i+1-i1;
            v = ae_c_mul(alpha,x->ptr.p_complex[i+1-i1]);
            ae_v_cmovec(&t->ptr.p_complex[tp1], 1, &y->ptr.p_complex[tp1], 1, "Conj", ae_v_len(tp1,tp2), v);
            v = ae_c_mul(ae_c_conj(alpha, _state),y->ptr.p_complex[i+1-i1]);
            ae_v_caddc(&t->ptr.p_complex[tp1], 1, &x->ptr.p_complex[tp1], 1, "Conj", ae_v_len(tp1,tp2), v);
            ae_v_cadd(&a->ptr.pp_complex[i][i1], 1, &t->ptr.p_complex[tp1], 1, "N", ae_v_len(i1,i));
        }
    }
}




/*************************************************************************
Generation of an elementary reflection transformation

The subroutine generates elementary reflection H of order N, so that, for
a given X, the following equality holds true:

    ( X(1) )   ( Beta )
H * (  ..  ) = (  0   )
    ( X(n) )   (  0   )

where
              ( V(1) )
H = 1 - Tau * (  ..  ) * ( V(1), ..., V(n) )
              ( V(n) )

where the first component of vector V equals 1.

Input parameters:
    X   -   vector. Array whose index ranges within [1..N].
    N   -   reflection order.

Output parameters:
    X   -   components from 2 to N are replaced with vector V.
            The first component is replaced with parameter Beta.
    Tau -   scalar value Tau. If X is a null vector, Tau equals 0,
            otherwise 1 <= Tau <= 2.

This subroutine is the modification of the DLARFG subroutines from
the LAPACK library.

MODIFICATIONS:
    24.12.2005 sign(Alpha) was replaced with an analogous to the Fortran SIGN code.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void generatereflection(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* tau,
     ae_state *_state)
{
    ae_int_t j;
    double alpha;
    double xnorm;
    double v;
    double beta;
    double mx;
    double s;

    *tau = 0;

    if( n<=1 )
    {
        *tau = (double)(0);
        return;
    }
    
    /*
     * Scale if needed (to avoid overflow/underflow during intermediate
     * calculations).
     */
    mx = (double)(0);
    for(j=1; j<=n; j++)
    {
        mx = ae_maxreal(ae_fabs(x->ptr.p_double[j], _state), mx, _state);
    }
    s = (double)(1);
    if( ae_fp_neq(mx,(double)(0)) )
    {
        if( ae_fp_less_eq(mx,ae_minrealnumber/ae_machineepsilon) )
        {
            s = ae_minrealnumber/ae_machineepsilon;
            v = 1/s;
            ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), v);
            mx = mx*v;
        }
        else
        {
            if( ae_fp_greater_eq(mx,ae_maxrealnumber*ae_machineepsilon) )
            {
                s = ae_maxrealnumber*ae_machineepsilon;
                v = 1/s;
                ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), v);
                mx = mx*v;
            }
        }
    }
    
    /*
     * XNORM = DNRM2( N-1, X, INCX )
     */
    alpha = x->ptr.p_double[1];
    xnorm = (double)(0);
    if( ae_fp_neq(mx,(double)(0)) )
    {
        for(j=2; j<=n; j++)
        {
            xnorm = xnorm+ae_sqr(x->ptr.p_double[j]/mx, _state);
        }
        xnorm = ae_sqrt(xnorm, _state)*mx;
    }
    if( ae_fp_eq(xnorm,(double)(0)) )
    {
        
        /*
         * H  =  I
         */
        *tau = (double)(0);
        x->ptr.p_double[1] = x->ptr.p_double[1]*s;
        return;
    }
    
    /*
     * general case
     */
    mx = ae_maxreal(ae_fabs(alpha, _state), ae_fabs(xnorm, _state), _state);
    beta = -mx*ae_sqrt(ae_sqr(alpha/mx, _state)+ae_sqr(xnorm/mx, _state), _state);
    if( ae_fp_less(alpha,(double)(0)) )
    {
        beta = -beta;
    }
    *tau = (beta-alpha)/beta;
    v = 1/(alpha-beta);
    ae_v_muld(&x->ptr.p_double[2], 1, ae_v_len(2,n), v);
    x->ptr.p_double[1] = beta;
    
    /*
     * Scale back outputs
     */
    x->ptr.p_double[1] = x->ptr.p_double[1]*s;
}


/*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The algorithm pre-multiplies the matrix by an elementary reflection transformation
which is given by column V and scalar Tau (see the description of the
GenerateReflection procedure). Not the whole matrix but only a part of it
is transformed (rows from M1 to M2, columns from N1 to N2). Only the elements
of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining the transformation.
    V       -   column defining the transformation.
                Array whose index ranges within [1..M2-M1+1].
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose indexes goes from N1 to N2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void applyreflectionfromtheleft(/* Real    */ ae_matrix* c,
     double tau,
     /* Real    */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* work,
     ae_state *_state)
{
    double t;
    ae_int_t i;


    if( (ae_fp_eq(tau,(double)(0))||n1>n2)||m1>m2 )
    {
        return;
    }
    
    /*
     * w := C' * v
     */
    for(i=n1; i<=n2; i++)
    {
        work->ptr.p_double[i] = (double)(0);
    }
    for(i=m1; i<=m2; i++)
    {
        t = v->ptr.p_double[i+1-m1];
        ae_v_addd(&work->ptr.p_double[n1], 1, &c->ptr.pp_double[i][n1], 1, ae_v_len(n1,n2), t);
    }
    
    /*
     * C := C - tau * v * w'
     */
    for(i=m1; i<=m2; i++)
    {
        t = v->ptr.p_double[i-m1+1]*tau;
        ae_v_subd(&c->ptr.pp_double[i][n1], 1, &work->ptr.p_double[n1], 1, ae_v_len(n1,n2), t);
    }
}


/*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The algorithm post-multiplies the matrix by an elementary reflection transformation
which is given by column V and scalar Tau (see the description of the
GenerateReflection procedure). Not the whole matrix but only a part of it
is transformed (rows from M1 to M2, columns from N1 to N2). Only the
elements of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining the transformation.
    V       -   column defining the transformation.
                Array whose index ranges within [1..N2-N1+1].
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose indexes goes from M1 to M2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void applyreflectionfromtheright(/* Real    */ ae_matrix* c,
     double tau,
     /* Real    */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* work,
     ae_state *_state)
{
    double t;
    ae_int_t i;
    ae_int_t vm;


    if( (ae_fp_eq(tau,(double)(0))||n1>n2)||m1>m2 )
    {
        return;
    }
    vm = n2-n1+1;
    for(i=m1; i<=m2; i++)
    {
        t = ae_v_dotproduct(&c->ptr.pp_double[i][n1], 1, &v->ptr.p_double[1], 1, ae_v_len(n1,n2));
        t = t*tau;
        ae_v_subd(&c->ptr.pp_double[i][n1], 1, &v->ptr.p_double[1], 1, ae_v_len(n1,n2), t);
    }
    
    /*
     * This line is necessary to avoid spurious compiler warnings
     */
    touchint(&vm, _state);
}




/*************************************************************************
Generation of an elementary complex reflection transformation

The subroutine generates elementary complex reflection H of  order  N,  so
that, for a given X, the following equality holds true:

     ( X(1) )   ( Beta )
H' * (  ..  ) = (  0   ),   H'*H = I,   Beta is a real number
     ( X(n) )   (  0   )

where

              ( V(1) )
H = 1 - Tau * (  ..  ) * ( conj(V(1)), ..., conj(V(n)) )
              ( V(n) )

where the first component of vector V equals 1.

Input parameters:
    X   -   vector. Array with elements [1..N].
    N   -   reflection order.

Output parameters:
    X   -   components from 2 to N are replaced by vector V.
            The first component is replaced with parameter Beta.
    Tau -   scalar value Tau.

This subroutine is the modification of CLARFG subroutines  from the LAPACK
library. It has similar functionality except for the fact that it  doesn’t
handle errors when intermediate results cause an overflow.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void complexgeneratereflection(/* Complex */ ae_vector* x,
     ae_int_t n,
     ae_complex* tau,
     ae_state *_state)
{
    ae_int_t j;
    ae_complex alpha;
    double alphi;
    double alphr;
    double beta;
    double xnorm;
    double mx;
    ae_complex t;
    double s;
    ae_complex v;

    tau->x = 0;
    tau->y = 0;

    if( n<=0 )
    {
        *tau = ae_complex_from_i(0);
        return;
    }
    
    /*
     * Scale if needed (to avoid overflow/underflow during intermediate
     * calculations).
     */
    mx = (double)(0);
    for(j=1; j<=n; j++)
    {
        mx = ae_maxreal(ae_c_abs(x->ptr.p_complex[j], _state), mx, _state);
    }
    s = (double)(1);
    if( ae_fp_neq(mx,(double)(0)) )
    {
        if( ae_fp_less(mx,(double)(1)) )
        {
            s = ae_sqrt(ae_minrealnumber, _state);
            v = ae_complex_from_d(1/s);
            ae_v_cmulc(&x->ptr.p_complex[1], 1, ae_v_len(1,n), v);
        }
        else
        {
            s = ae_sqrt(ae_maxrealnumber, _state);
            v = ae_complex_from_d(1/s);
            ae_v_cmulc(&x->ptr.p_complex[1], 1, ae_v_len(1,n), v);
        }
    }
    
    /*
     * calculate
     */
    alpha = x->ptr.p_complex[1];
    mx = (double)(0);
    for(j=2; j<=n; j++)
    {
        mx = ae_maxreal(ae_c_abs(x->ptr.p_complex[j], _state), mx, _state);
    }
    xnorm = (double)(0);
    if( ae_fp_neq(mx,(double)(0)) )
    {
        for(j=2; j<=n; j++)
        {
            t = ae_c_div_d(x->ptr.p_complex[j],mx);
            xnorm = xnorm+ae_c_mul(t,ae_c_conj(t, _state)).x;
        }
        xnorm = ae_sqrt(xnorm, _state)*mx;
    }
    alphr = alpha.x;
    alphi = alpha.y;
    if( ae_fp_eq(xnorm,(double)(0))&&ae_fp_eq(alphi,(double)(0)) )
    {
        *tau = ae_complex_from_i(0);
        x->ptr.p_complex[1] = ae_c_mul_d(x->ptr.p_complex[1],s);
        return;
    }
    mx = ae_maxreal(ae_fabs(alphr, _state), ae_fabs(alphi, _state), _state);
    mx = ae_maxreal(mx, ae_fabs(xnorm, _state), _state);
    beta = -mx*ae_sqrt(ae_sqr(alphr/mx, _state)+ae_sqr(alphi/mx, _state)+ae_sqr(xnorm/mx, _state), _state);
    if( ae_fp_less(alphr,(double)(0)) )
    {
        beta = -beta;
    }
    tau->x = (beta-alphr)/beta;
    tau->y = -alphi/beta;
    alpha = ae_c_d_div(1,ae_c_sub_d(alpha,beta));
    if( n>1 )
    {
        ae_v_cmulc(&x->ptr.p_complex[2], 1, ae_v_len(2,n), alpha);
    }
    alpha = ae_complex_from_d(beta);
    x->ptr.p_complex[1] = alpha;
    
    /*
     * Scale back
     */
    x->ptr.p_complex[1] = ae_c_mul_d(x->ptr.p_complex[1],s);
}


/*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The  algorithm  pre-multiplies  the  matrix  by  an  elementary reflection
transformation  which  is  given  by  column  V  and  scalar  Tau (see the
description of the GenerateReflection). Not the whole matrix  but  only  a
part of it is transformed (rows from M1 to M2, columns from N1 to N2). Only
the elements of this submatrix are changed.

Note: the matrix is multiplied by H, not by H'.   If  it  is  required  to
multiply the matrix by H', it is necessary to pass Conj(Tau) instead of Tau.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining transformation.
    V       -   column defining transformation.
                Array whose index ranges within [1..M2-M1+1]
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose index goes from N1 to N2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void complexapplyreflectionfromtheleft(/* Complex */ ae_matrix* c,
     ae_complex tau,
     /* Complex */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Complex */ ae_vector* work,
     ae_state *_state)
{
    ae_complex t;
    ae_int_t i;


    if( (ae_c_eq_d(tau,(double)(0))||n1>n2)||m1>m2 )
    {
        return;
    }
    
    /*
     * w := C^T * conj(v)
     */
    for(i=n1; i<=n2; i++)
    {
        work->ptr.p_complex[i] = ae_complex_from_i(0);
    }
    for(i=m1; i<=m2; i++)
    {
        t = ae_c_conj(v->ptr.p_complex[i+1-m1], _state);
        ae_v_caddc(&work->ptr.p_complex[n1], 1, &c->ptr.pp_complex[i][n1], 1, "N", ae_v_len(n1,n2), t);
    }
    
    /*
     * C := C - tau * v * w^T
     */
    for(i=m1; i<=m2; i++)
    {
        t = ae_c_mul(v->ptr.p_complex[i-m1+1],tau);
        ae_v_csubc(&c->ptr.pp_complex[i][n1], 1, &work->ptr.p_complex[n1], 1, "N", ae_v_len(n1,n2), t);
    }
}


/*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The  algorithm  post-multiplies  the  matrix  by  an elementary reflection
transformation  which  is  given  by  column  V  and  scalar  Tau (see the
description  of  the  GenerateReflection). Not the whole matrix but only a
part  of  it  is  transformed (rows from M1 to M2, columns from N1 to N2).
Only the elements of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining transformation.
    V       -   column defining transformation.
                Array whose index ranges within [1..N2-N1+1]
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose index goes from M1 to M2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void complexapplyreflectionfromtheright(/* Complex */ ae_matrix* c,
     ae_complex tau,
     /* Complex */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Complex */ ae_vector* work,
     ae_state *_state)
{
    ae_complex t;
    ae_int_t i;
    ae_int_t vm;


    if( (ae_c_eq_d(tau,(double)(0))||n1>n2)||m1>m2 )
    {
        return;
    }
    
    /*
     * w := C * v
     */
    vm = n2-n1+1;
    for(i=m1; i<=m2; i++)
    {
        t = ae_v_cdotproduct(&c->ptr.pp_complex[i][n1], 1, "N", &v->ptr.p_complex[1], 1, "N", ae_v_len(n1,n2));
        work->ptr.p_complex[i] = t;
    }
    
    /*
     * C := C - w * conj(v^T)
     */
    ae_v_cmove(&v->ptr.p_complex[1], 1, &v->ptr.p_complex[1], 1, "Conj", ae_v_len(1,vm));
    for(i=m1; i<=m2; i++)
    {
        t = ae_c_mul(work->ptr.p_complex[i],tau);
        ae_v_csubc(&c->ptr.pp_complex[i][n1], 1, &v->ptr.p_complex[1], 1, "N", ae_v_len(n1,n2), t);
    }
    ae_v_cmove(&v->ptr.p_complex[1], 1, &v->ptr.p_complex[1], 1, "Conj", ae_v_len(1,vm));
}




void symmetricmatrixvectormultiply(/* Real    */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Real    */ ae_vector* x,
     double alpha,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ba1;
    ae_int_t ba2;
    ae_int_t by1;
    ae_int_t by2;
    ae_int_t bx1;
    ae_int_t bx2;
    ae_int_t n;
    double v;


    n = i2-i1+1;
    if( n<=0 )
    {
        return;
    }
    
    /*
     * Let A = L + D + U, where
     *  L is strictly lower triangular (main diagonal is zero)
     *  D is diagonal
     *  U is strictly upper triangular (main diagonal is zero)
     *
     * A*x = L*x + D*x + U*x
     *
     * Calculate D*x first
     */
    for(i=i1; i<=i2; i++)
    {
        y->ptr.p_double[i-i1+1] = a->ptr.pp_double[i][i]*x->ptr.p_double[i-i1+1];
    }
    
    /*
     * Add L*x + U*x
     */
    if( isupper )
    {
        for(i=i1; i<=i2-1; i++)
        {
            
            /*
             * Add L*x to the result
             */
            v = x->ptr.p_double[i-i1+1];
            by1 = i-i1+2;
            by2 = n;
            ba1 = i+1;
            ba2 = i2;
            ae_v_addd(&y->ptr.p_double[by1], 1, &a->ptr.pp_double[i][ba1], 1, ae_v_len(by1,by2), v);
            
            /*
             * Add U*x to the result
             */
            bx1 = i-i1+2;
            bx2 = n;
            ba1 = i+1;
            ba2 = i2;
            v = ae_v_dotproduct(&x->ptr.p_double[bx1], 1, &a->ptr.pp_double[i][ba1], 1, ae_v_len(bx1,bx2));
            y->ptr.p_double[i-i1+1] = y->ptr.p_double[i-i1+1]+v;
        }
    }
    else
    {
        for(i=i1+1; i<=i2; i++)
        {
            
            /*
             * Add L*x to the result
             */
            bx1 = 1;
            bx2 = i-i1;
            ba1 = i1;
            ba2 = i-1;
            v = ae_v_dotproduct(&x->ptr.p_double[bx1], 1, &a->ptr.pp_double[i][ba1], 1, ae_v_len(bx1,bx2));
            y->ptr.p_double[i-i1+1] = y->ptr.p_double[i-i1+1]+v;
            
            /*
             * Add U*x to the result
             */
            v = x->ptr.p_double[i-i1+1];
            by1 = 1;
            by2 = i-i1;
            ba1 = i1;
            ba2 = i-1;
            ae_v_addd(&y->ptr.p_double[by1], 1, &a->ptr.pp_double[i][ba1], 1, ae_v_len(by1,by2), v);
        }
    }
    ae_v_muld(&y->ptr.p_double[1], 1, ae_v_len(1,n), alpha);
    touchint(&ba2, _state);
}


void symmetricrank2update(/* Real    */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* t,
     double alpha,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t tp1;
    ae_int_t tp2;
    double v;


    if( isupper )
    {
        for(i=i1; i<=i2; i++)
        {
            tp1 = i+1-i1;
            tp2 = i2-i1+1;
            v = x->ptr.p_double[i+1-i1];
            ae_v_moved(&t->ptr.p_double[tp1], 1, &y->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), v);
            v = y->ptr.p_double[i+1-i1];
            ae_v_addd(&t->ptr.p_double[tp1], 1, &x->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), v);
            ae_v_muld(&t->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), alpha);
            ae_v_add(&a->ptr.pp_double[i][i], 1, &t->ptr.p_double[tp1], 1, ae_v_len(i,i2));
        }
    }
    else
    {
        for(i=i1; i<=i2; i++)
        {
            tp1 = 1;
            tp2 = i+1-i1;
            v = x->ptr.p_double[i+1-i1];
            ae_v_moved(&t->ptr.p_double[tp1], 1, &y->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), v);
            v = y->ptr.p_double[i+1-i1];
            ae_v_addd(&t->ptr.p_double[tp1], 1, &x->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), v);
            ae_v_muld(&t->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), alpha);
            ae_v_add(&a->ptr.pp_double[i][i1], 1, &t->ptr.p_double[tp1], 1, ae_v_len(i1,i));
        }
    }
}




/*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm pre-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on, are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1 to
M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..M2-M1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [N1..N2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************/
void applyrotationsfromtheleft(ae_bool isforward,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* c,
     /* Real    */ ae_vector* s,
     /* Real    */ ae_matrix* a,
     /* Real    */ ae_vector* work,
     ae_state *_state)
{
    ae_int_t j;
    ae_int_t jp1;
    double ctemp;
    double stemp;
    double temp;


    if( m1>m2||n1>n2 )
    {
        return;
    }
    
    /*
     * Form  P * A
     */
    if( isforward )
    {
        if( n1!=n2 )
        {
            
            /*
             * Common case: N1<>N2
             */
            for(j=m1; j<=m2-1; j++)
            {
                ctemp = c->ptr.p_double[j-m1+1];
                stemp = s->ptr.p_double[j-m1+1];
                if( ae_fp_neq(ctemp,(double)(1))||ae_fp_neq(stemp,(double)(0)) )
                {
                    jp1 = j+1;
                    ae_v_moved(&work->ptr.p_double[n1], 1, &a->ptr.pp_double[jp1][n1], 1, ae_v_len(n1,n2), ctemp);
                    ae_v_subd(&work->ptr.p_double[n1], 1, &a->ptr.pp_double[j][n1], 1, ae_v_len(n1,n2), stemp);
                    ae_v_muld(&a->ptr.pp_double[j][n1], 1, ae_v_len(n1,n2), ctemp);
                    ae_v_addd(&a->ptr.pp_double[j][n1], 1, &a->ptr.pp_double[jp1][n1], 1, ae_v_len(n1,n2), stemp);
                    ae_v_move(&a->ptr.pp_double[jp1][n1], 1, &work->ptr.p_double[n1], 1, ae_v_len(n1,n2));
                }
            }
        }
        else
        {
            
            /*
             * Special case: N1=N2
             */
            for(j=m1; j<=m2-1; j++)
            {
                ctemp = c->ptr.p_double[j-m1+1];
                stemp = s->ptr.p_double[j-m1+1];
                if( ae_fp_neq(ctemp,(double)(1))||ae_fp_neq(stemp,(double)(0)) )
                {
                    temp = a->ptr.pp_double[j+1][n1];
                    a->ptr.pp_double[j+1][n1] = ctemp*temp-stemp*a->ptr.pp_double[j][n1];
                    a->ptr.pp_double[j][n1] = stemp*temp+ctemp*a->ptr.pp_double[j][n1];
                }
            }
        }
    }
    else
    {
        if( n1!=n2 )
        {
            
            /*
             * Common case: N1<>N2
             */
            for(j=m2-1; j>=m1; j--)
            {
                ctemp = c->ptr.p_double[j-m1+1];
                stemp = s->ptr.p_double[j-m1+1];
                if( ae_fp_neq(ctemp,(double)(1))||ae_fp_neq(stemp,(double)(0)) )
                {
                    jp1 = j+1;
                    ae_v_moved(&work->ptr.p_double[n1], 1, &a->ptr.pp_double[jp1][n1], 1, ae_v_len(n1,n2), ctemp);
                    ae_v_subd(&work->ptr.p_double[n1], 1, &a->ptr.pp_double[j][n1], 1, ae_v_len(n1,n2), stemp);
                    ae_v_muld(&a->ptr.pp_double[j][n1], 1, ae_v_len(n1,n2), ctemp);
                    ae_v_addd(&a->ptr.pp_double[j][n1], 1, &a->ptr.pp_double[jp1][n1], 1, ae_v_len(n1,n2), stemp);
                    ae_v_move(&a->ptr.pp_double[jp1][n1], 1, &work->ptr.p_double[n1], 1, ae_v_len(n1,n2));
                }
            }
        }
        else
        {
            
            /*
             * Special case: N1=N2
             */
            for(j=m2-1; j>=m1; j--)
            {
                ctemp = c->ptr.p_double[j-m1+1];
                stemp = s->ptr.p_double[j-m1+1];
                if( ae_fp_neq(ctemp,(double)(1))||ae_fp_neq(stemp,(double)(0)) )
                {
                    temp = a->ptr.pp_double[j+1][n1];
                    a->ptr.pp_double[j+1][n1] = ctemp*temp-stemp*a->ptr.pp_double[j][n1];
                    a->ptr.pp_double[j][n1] = stemp*temp+ctemp*a->ptr.pp_double[j][n1];
                }
            }
        }
    }
}


/*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm post-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1
to M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..N2-N1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [M1..M2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************/
void applyrotationsfromtheright(ae_bool isforward,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* c,
     /* Real    */ ae_vector* s,
     /* Real    */ ae_matrix* a,
     /* Real    */ ae_vector* work,
     ae_state *_state)
{
    ae_int_t j;
    ae_int_t jp1;
    double ctemp;
    double stemp;
    double temp;


    
    /*
     * Form A * P'
     */
    if( isforward )
    {
        if( m1!=m2 )
        {
            
            /*
             * Common case: M1<>M2
             */
            for(j=n1; j<=n2-1; j++)
            {
                ctemp = c->ptr.p_double[j-n1+1];
                stemp = s->ptr.p_double[j-n1+1];
                if( ae_fp_neq(ctemp,(double)(1))||ae_fp_neq(stemp,(double)(0)) )
                {
                    jp1 = j+1;
                    ae_v_moved(&work->ptr.p_double[m1], 1, &a->ptr.pp_double[m1][jp1], a->stride, ae_v_len(m1,m2), ctemp);
                    ae_v_subd(&work->ptr.p_double[m1], 1, &a->ptr.pp_double[m1][j], a->stride, ae_v_len(m1,m2), stemp);
                    ae_v_muld(&a->ptr.pp_double[m1][j], a->stride, ae_v_len(m1,m2), ctemp);
                    ae_v_addd(&a->ptr.pp_double[m1][j], a->stride, &a->ptr.pp_double[m1][jp1], a->stride, ae_v_len(m1,m2), stemp);
                    ae_v_move(&a->ptr.pp_double[m1][jp1], a->stride, &work->ptr.p_double[m1], 1, ae_v_len(m1,m2));
                }
            }
        }
        else
        {
            
            /*
             * Special case: M1=M2
             */
            for(j=n1; j<=n2-1; j++)
            {
                ctemp = c->ptr.p_double[j-n1+1];
                stemp = s->ptr.p_double[j-n1+1];
                if( ae_fp_neq(ctemp,(double)(1))||ae_fp_neq(stemp,(double)(0)) )
                {
                    temp = a->ptr.pp_double[m1][j+1];
                    a->ptr.pp_double[m1][j+1] = ctemp*temp-stemp*a->ptr.pp_double[m1][j];
                    a->ptr.pp_double[m1][j] = stemp*temp+ctemp*a->ptr.pp_double[m1][j];
                }
            }
        }
    }
    else
    {
        if( m1!=m2 )
        {
            
            /*
             * Common case: M1<>M2
             */
            for(j=n2-1; j>=n1; j--)
            {
                ctemp = c->ptr.p_double[j-n1+1];
                stemp = s->ptr.p_double[j-n1+1];
                if( ae_fp_neq(ctemp,(double)(1))||ae_fp_neq(stemp,(double)(0)) )
                {
                    jp1 = j+1;
                    ae_v_moved(&work->ptr.p_double[m1], 1, &a->ptr.pp_double[m1][jp1], a->stride, ae_v_len(m1,m2), ctemp);
                    ae_v_subd(&work->ptr.p_double[m1], 1, &a->ptr.pp_double[m1][j], a->stride, ae_v_len(m1,m2), stemp);
                    ae_v_muld(&a->ptr.pp_double[m1][j], a->stride, ae_v_len(m1,m2), ctemp);
                    ae_v_addd(&a->ptr.pp_double[m1][j], a->stride, &a->ptr.pp_double[m1][jp1], a->stride, ae_v_len(m1,m2), stemp);
                    ae_v_move(&a->ptr.pp_double[m1][jp1], a->stride, &work->ptr.p_double[m1], 1, ae_v_len(m1,m2));
                }
            }
        }
        else
        {
            
            /*
             * Special case: M1=M2
             */
            for(j=n2-1; j>=n1; j--)
            {
                ctemp = c->ptr.p_double[j-n1+1];
                stemp = s->ptr.p_double[j-n1+1];
                if( ae_fp_neq(ctemp,(double)(1))||ae_fp_neq(stemp,(double)(0)) )
                {
                    temp = a->ptr.pp_double[m1][j+1];
                    a->ptr.pp_double[m1][j+1] = ctemp*temp-stemp*a->ptr.pp_double[m1][j];
                    a->ptr.pp_double[m1][j] = stemp*temp+ctemp*a->ptr.pp_double[m1][j];
                }
            }
        }
    }
}


/*************************************************************************
The subroutine generates the elementary rotation, so that:

[  CS  SN  ]  .  [ F ]  =  [ R ]
[ -SN  CS  ]     [ G ]     [ 0 ]

CS**2 + SN**2 = 1
*************************************************************************/
void generaterotation(double f,
     double g,
     double* cs,
     double* sn,
     double* r,
     ae_state *_state)
{
    double f1;
    double g1;

    *cs = 0;
    *sn = 0;
    *r = 0;

    if( ae_fp_eq(g,(double)(0)) )
    {
        *cs = (double)(1);
        *sn = (double)(0);
        *r = f;
    }
    else
    {
        if( ae_fp_eq(f,(double)(0)) )
        {
            *cs = (double)(0);
            *sn = (double)(1);
            *r = g;
        }
        else
        {
            f1 = f;
            g1 = g;
            if( ae_fp_greater(ae_fabs(f1, _state),ae_fabs(g1, _state)) )
            {
                *r = ae_fabs(f1, _state)*ae_sqrt(1+ae_sqr(g1/f1, _state), _state);
            }
            else
            {
                *r = ae_fabs(g1, _state)*ae_sqrt(1+ae_sqr(f1/g1, _state), _state);
            }
            *cs = f1/(*r);
            *sn = g1/(*r);
            if( ae_fp_greater(ae_fabs(f, _state),ae_fabs(g, _state))&&ae_fp_less(*cs,(double)(0)) )
            {
                *cs = -*cs;
                *sn = -*sn;
                *r = -*r;
            }
        }
    }
}




void rmatrixinternalschurdecomposition(/* Real    */ ae_matrix* h,
     ae_int_t n,
     ae_int_t tneeded,
     ae_int_t zneeded,
     /* Real    */ ae_vector* wr,
     /* Real    */ ae_vector* wi,
     /* Real    */ ae_matrix* z,
     ae_int_t* info,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_matrix h1;
    ae_matrix z1;
    ae_vector wr1;
    ae_vector wi1;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(wr);
    ae_vector_clear(wi);
    *info = 0;
    ae_matrix_init(&h1, 0, 0, DT_REAL, _state);
    ae_matrix_init(&z1, 0, 0, DT_REAL, _state);
    ae_vector_init(&wr1, 0, DT_REAL, _state);
    ae_vector_init(&wi1, 0, DT_REAL, _state);

    
    /*
     * Allocate space
     */
    ae_vector_set_length(wr, n, _state);
    ae_vector_set_length(wi, n, _state);
    if( zneeded==2 )
    {
        rmatrixsetlengthatleast(z, n, n, _state);
    }
    
    /*
     * MKL version
     */
    if( rmatrixinternalschurdecompositionmkl(h, n, tneeded, zneeded, wr, wi, z, info, _state) )
    {
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * ALGLIB version
     */
    ae_matrix_set_length(&h1, n+1, n+1, _state);
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            h1.ptr.pp_double[1+i][1+j] = h->ptr.pp_double[i][j];
        }
    }
    if( zneeded==1 )
    {
        ae_matrix_set_length(&z1, n+1, n+1, _state);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                z1.ptr.pp_double[1+i][1+j] = z->ptr.pp_double[i][j];
            }
        }
    }
    internalschurdecomposition(&h1, n, tneeded, zneeded, &wr1, &wi1, &z1, info, _state);
    for(i=0; i<=n-1; i++)
    {
        wr->ptr.p_double[i] = wr1.ptr.p_double[i+1];
        wi->ptr.p_double[i] = wi1.ptr.p_double[i+1];
    }
    if( tneeded!=0 )
    {
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                h->ptr.pp_double[i][j] = h1.ptr.pp_double[1+i][1+j];
            }
        }
    }
    if( zneeded!=0 )
    {
        rmatrixsetlengthatleast(z, n, n, _state);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                z->ptr.pp_double[i][j] = z1.ptr.pp_double[1+i][1+j];
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Subroutine performing  the  Schur  decomposition  of  a  matrix  in  upper
Hessenberg form using the QR algorithm with multiple shifts.

The  source matrix  H  is  represented as  S'*H*S = T, where H - matrix in
upper Hessenberg form,  S - orthogonal matrix (Schur vectors),   T - upper
quasi-triangular matrix (with blocks of sizes  1x1  and  2x2  on  the main
diagonal).

Input parameters:
    H   -   matrix to be decomposed.
            Array whose indexes range within [1..N, 1..N].
    N   -   size of H, N>=0.


Output parameters:
    H   –   contains the matrix T.
            Array whose indexes range within [1..N, 1..N].
            All elements below the blocks on the main diagonal are equal
            to 0.
    S   -   contains Schur vectors.
            Array whose indexes range within [1..N, 1..N].

Note 1:
    The block structure of matrix T could be easily recognized: since  all
    the elements  below  the blocks are zeros, the elements a[i+1,i] which
    are equal to 0 show the block border.

Note 2:
    the algorithm  performance  depends  on  the  value  of  the  internal
    parameter NS of InternalSchurDecomposition  subroutine  which  defines
    the number of shifts in the QR algorithm (analog of  the  block  width
    in block matrix algorithms in linear algebra). If you require  maximum
    performance  on  your  machine,  it  is  recommended  to  adjust  this
    parameter manually.

Result:
    True, if the algorithm has converged and the parameters H and S contain
        the result.
    False, if the algorithm has not converged.

Algorithm implemented on the basis of subroutine DHSEQR (LAPACK 3.0 library).
*************************************************************************/
ae_bool upperhessenbergschurdecomposition(/* Real    */ ae_matrix* h,
     ae_int_t n,
     /* Real    */ ae_matrix* s,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector wi;
    ae_vector wr;
    ae_int_t info;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_clear(s);
    ae_vector_init(&wi, 0, DT_REAL, _state);
    ae_vector_init(&wr, 0, DT_REAL, _state);

    internalschurdecomposition(h, n, 1, 2, &wr, &wi, s, &info, _state);
    result = info==0;
    ae_frame_leave(_state);
    return result;
}


void internalschurdecomposition(/* Real    */ ae_matrix* h,
     ae_int_t n,
     ae_int_t tneeded,
     ae_int_t zneeded,
     /* Real    */ ae_vector* wr,
     /* Real    */ ae_vector* wi,
     /* Real    */ ae_matrix* z,
     ae_int_t* info,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector work;
    ae_int_t i;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t ierr;
    ae_int_t ii;
    ae_int_t itemp;
    ae_int_t itn;
    ae_int_t its;
    ae_int_t j;
    ae_int_t k;
    ae_int_t l;
    ae_int_t maxb;
    ae_int_t nr;
    ae_int_t ns;
    ae_int_t nv;
    double absw;
    double smlnum;
    double tau;
    double temp;
    double tst1;
    double ulp;
    double unfl;
    ae_matrix s;
    ae_vector v;
    ae_vector vv;
    ae_vector workc1;
    ae_vector works1;
    ae_vector workv3;
    ae_vector tmpwr;
    ae_vector tmpwi;
    ae_bool initz;
    ae_bool wantt;
    ae_bool wantz;
    double cnst;
    ae_bool failflag;
    ae_int_t p1;
    ae_int_t p2;
    double vt;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(wr);
    ae_vector_clear(wi);
    *info = 0;
    ae_vector_init(&work, 0, DT_REAL, _state);
    ae_matrix_init(&s, 0, 0, DT_REAL, _state);
    ae_vector_init(&v, 0, DT_REAL, _state);
    ae_vector_init(&vv, 0, DT_REAL, _state);
    ae_vector_init(&workc1, 0, DT_REAL, _state);
    ae_vector_init(&works1, 0, DT_REAL, _state);
    ae_vector_init(&workv3, 0, DT_REAL, _state);
    ae_vector_init(&tmpwr, 0, DT_REAL, _state);
    ae_vector_init(&tmpwi, 0, DT_REAL, _state);

    
    /*
     * Set the order of the multi-shift QR algorithm to be used.
     * If you want to tune algorithm, change this values
     */
    ns = 12;
    maxb = 50;
    
    /*
     * Now 2 < NS <= MAXB < NH.
     */
    maxb = ae_maxint(3, maxb, _state);
    ns = ae_minint(maxb, ns, _state);
    
    /*
     * Initialize
     */
    cnst = 1.5;
    ae_vector_set_length(&work, ae_maxint(n, 1, _state)+1, _state);
    ae_matrix_set_length(&s, ns+1, ns+1, _state);
    ae_vector_set_length(&v, ns+1+1, _state);
    ae_vector_set_length(&vv, ns+1+1, _state);
    ae_vector_set_length(wr, ae_maxint(n, 1, _state)+1, _state);
    ae_vector_set_length(wi, ae_maxint(n, 1, _state)+1, _state);
    ae_vector_set_length(&workc1, 1+1, _state);
    ae_vector_set_length(&works1, 1+1, _state);
    ae_vector_set_length(&workv3, 3+1, _state);
    ae_vector_set_length(&tmpwr, ae_maxint(n, 1, _state)+1, _state);
    ae_vector_set_length(&tmpwi, ae_maxint(n, 1, _state)+1, _state);
    ae_assert(n>=0, "InternalSchurDecomposition: incorrect N!", _state);
    ae_assert(tneeded==0||tneeded==1, "InternalSchurDecomposition: incorrect TNeeded!", _state);
    ae_assert((zneeded==0||zneeded==1)||zneeded==2, "InternalSchurDecomposition: incorrect ZNeeded!", _state);
    wantt = tneeded==1;
    initz = zneeded==2;
    wantz = zneeded!=0;
    *info = 0;
    
    /*
     * Initialize Z, if necessary
     */
    if( initz )
    {
        rmatrixsetlengthatleast(z, n+1, n+1, _state);
        for(i=1; i<=n; i++)
        {
            for(j=1; j<=n; j++)
            {
                if( i==j )
                {
                    z->ptr.pp_double[i][j] = (double)(1);
                }
                else
                {
                    z->ptr.pp_double[i][j] = (double)(0);
                }
            }
        }
    }
    
    /*
     * Quick return if possible
     */
    if( n==0 )
    {
        ae_frame_leave(_state);
        return;
    }
    if( n==1 )
    {
        wr->ptr.p_double[1] = h->ptr.pp_double[1][1];
        wi->ptr.p_double[1] = (double)(0);
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Set rows and columns 1 to N to zero below the first
     * subdiagonal.
     */
    for(j=1; j<=n-2; j++)
    {
        for(i=j+2; i<=n; i++)
        {
            h->ptr.pp_double[i][j] = (double)(0);
        }
    }
    
    /*
     * Test if N is sufficiently small
     */
    if( (ns<=2||ns>n)||maxb>=n )
    {
        
        /*
         * Use the standard double-shift algorithm
         */
        hsschur_internalauxschur(wantt, wantz, n, 1, n, h, wr, wi, 1, n, z, &work, &workv3, &workc1, &works1, info, _state);
        
        /*
         * fill entries under diagonal blocks of T with zeros
         */
        if( wantt )
        {
            j = 1;
            while(j<=n)
            {
                if( ae_fp_eq(wi->ptr.p_double[j],(double)(0)) )
                {
                    for(i=j+1; i<=n; i++)
                    {
                        h->ptr.pp_double[i][j] = (double)(0);
                    }
                    j = j+1;
                }
                else
                {
                    for(i=j+2; i<=n; i++)
                    {
                        h->ptr.pp_double[i][j] = (double)(0);
                        h->ptr.pp_double[i][j+1] = (double)(0);
                    }
                    j = j+2;
                }
            }
        }
        ae_frame_leave(_state);
        return;
    }
    unfl = ae_minrealnumber;
    ulp = 2*ae_machineepsilon;
    smlnum = unfl*(n/ulp);
    
    /*
     * I1 and I2 are the indices of the first row and last column of H
     * to which transformations must be applied. If eigenvalues only are
     * being computed, I1 and I2 are set inside the main loop.
     */
    i1 = 1;
    i2 = n;
    
    /*
     * ITN is the total number of multiple-shift QR iterations allowed.
     */
    itn = 30*n;
    
    /*
     * The main loop begins here. I is the loop index and decreases from
     * IHI to ILO in steps of at most MAXB. Each iteration of the loop
     * works with the active submatrix in rows and columns L to I.
     * Eigenvalues I+1 to IHI have already converged. Either L = ILO or
     * H(L,L-1) is negligible so that the matrix splits.
     */
    i = n;
    for(;;)
    {
        l = 1;
        if( i<1 )
        {
            
            /*
             * fill entries under diagonal blocks of T with zeros
             */
            if( wantt )
            {
                j = 1;
                while(j<=n)
                {
                    if( ae_fp_eq(wi->ptr.p_double[j],(double)(0)) )
                    {
                        for(i=j+1; i<=n; i++)
                        {
                            h->ptr.pp_double[i][j] = (double)(0);
                        }
                        j = j+1;
                    }
                    else
                    {
                        for(i=j+2; i<=n; i++)
                        {
                            h->ptr.pp_double[i][j] = (double)(0);
                            h->ptr.pp_double[i][j+1] = (double)(0);
                        }
                        j = j+2;
                    }
                }
            }
            
            /*
             * Exit
             */
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * Perform multiple-shift QR iterations on rows and columns ILO to I
         * until a submatrix of order at most MAXB splits off at the bottom
         * because a subdiagonal element has become negligible.
         */
        failflag = ae_true;
        for(its=0; its<=itn; its++)
        {
            
            /*
             * Look for a single small subdiagonal element.
             */
            for(k=i; k>=l+1; k--)
            {
                tst1 = ae_fabs(h->ptr.pp_double[k-1][k-1], _state)+ae_fabs(h->ptr.pp_double[k][k], _state);
                if( ae_fp_eq(tst1,(double)(0)) )
                {
                    tst1 = upperhessenberg1norm(h, l, i, l, i, &work, _state);
                }
                if( ae_fp_less_eq(ae_fabs(h->ptr.pp_double[k][k-1], _state),ae_maxreal(ulp*tst1, smlnum, _state)) )
                {
                    break;
                }
            }
            l = k;
            if( l>1 )
            {
                
                /*
                 * H(L,L-1) is negligible.
                 */
                h->ptr.pp_double[l][l-1] = (double)(0);
            }
            
            /*
             * Exit from loop if a submatrix of order <= MAXB has split off.
             */
            if( l>=i-maxb+1 )
            {
                failflag = ae_false;
                break;
            }
            
            /*
             * Now the active submatrix is in rows and columns L to I. If
             * eigenvalues only are being computed, only the active submatrix
             * need be transformed.
             */
            if( its==20||its==30 )
            {
                
                /*
                 * Exceptional shifts.
                 */
                for(ii=i-ns+1; ii<=i; ii++)
                {
                    wr->ptr.p_double[ii] = cnst*(ae_fabs(h->ptr.pp_double[ii][ii-1], _state)+ae_fabs(h->ptr.pp_double[ii][ii], _state));
                    wi->ptr.p_double[ii] = (double)(0);
                }
            }
            else
            {
                
                /*
                 * Use eigenvalues of trailing submatrix of order NS as shifts.
                 */
                copymatrix(h, i-ns+1, i, i-ns+1, i, &s, 1, ns, 1, ns, _state);
                hsschur_internalauxschur(ae_false, ae_false, ns, 1, ns, &s, &tmpwr, &tmpwi, 1, ns, z, &work, &workv3, &workc1, &works1, &ierr, _state);
                for(p1=1; p1<=ns; p1++)
                {
                    wr->ptr.p_double[i-ns+p1] = tmpwr.ptr.p_double[p1];
                    wi->ptr.p_double[i-ns+p1] = tmpwi.ptr.p_double[p1];
                }
                if( ierr>0 )
                {
                    
                    /*
                     * If DLAHQR failed to compute all NS eigenvalues, use the
                     * unconverged diagonal elements as the remaining shifts.
                     */
                    for(ii=1; ii<=ierr; ii++)
                    {
                        wr->ptr.p_double[i-ns+ii] = s.ptr.pp_double[ii][ii];
                        wi->ptr.p_double[i-ns+ii] = (double)(0);
                    }
                }
            }
            
            /*
             * Form the first column of (G-w(1)) (G-w(2)) . . . (G-w(ns))
             * where G is the Hessenberg submatrix H(L:I,L:I) and w is
             * the vector of shifts (stored in WR and WI). The result is
             * stored in the local array V.
             */
            v.ptr.p_double[1] = (double)(1);
            for(ii=2; ii<=ns+1; ii++)
            {
                v.ptr.p_double[ii] = (double)(0);
            }
            nv = 1;
            for(j=i-ns+1; j<=i; j++)
            {
                if( ae_fp_greater_eq(wi->ptr.p_double[j],(double)(0)) )
                {
                    if( ae_fp_eq(wi->ptr.p_double[j],(double)(0)) )
                    {
                        
                        /*
                         * real shift
                         */
                        p1 = nv+1;
                        ae_v_move(&vv.ptr.p_double[1], 1, &v.ptr.p_double[1], 1, ae_v_len(1,p1));
                        matrixvectormultiply(h, l, l+nv, l, l+nv-1, ae_false, &vv, 1, nv, 1.0, &v, 1, nv+1, -wr->ptr.p_double[j], _state);
                        nv = nv+1;
                    }
                    else
                    {
                        if( ae_fp_greater(wi->ptr.p_double[j],(double)(0)) )
                        {
                            
                            /*
                             * complex conjugate pair of shifts
                             */
                            p1 = nv+1;
                            ae_v_move(&vv.ptr.p_double[1], 1, &v.ptr.p_double[1], 1, ae_v_len(1,p1));
                            matrixvectormultiply(h, l, l+nv, l, l+nv-1, ae_false, &v, 1, nv, 1.0, &vv, 1, nv+1, -2*wr->ptr.p_double[j], _state);
                            itemp = vectoridxabsmax(&vv, 1, nv+1, _state);
                            temp = 1/ae_maxreal(ae_fabs(vv.ptr.p_double[itemp], _state), smlnum, _state);
                            p1 = nv+1;
                            ae_v_muld(&vv.ptr.p_double[1], 1, ae_v_len(1,p1), temp);
                            absw = pythag2(wr->ptr.p_double[j], wi->ptr.p_double[j], _state);
                            temp = temp*absw*absw;
                            matrixvectormultiply(h, l, l+nv+1, l, l+nv, ae_false, &vv, 1, nv+1, 1.0, &v, 1, nv+2, temp, _state);
                            nv = nv+2;
                        }
                    }
                    
                    /*
                     * Scale V(1:NV) so that max(abs(V(i))) = 1. If V is zero,
                     * reset it to the unit vector.
                     */
                    itemp = vectoridxabsmax(&v, 1, nv, _state);
                    temp = ae_fabs(v.ptr.p_double[itemp], _state);
                    if( ae_fp_eq(temp,(double)(0)) )
                    {
                        v.ptr.p_double[1] = (double)(1);
                        for(ii=2; ii<=nv; ii++)
                        {
                            v.ptr.p_double[ii] = (double)(0);
                        }
                    }
                    else
                    {
                        temp = ae_maxreal(temp, smlnum, _state);
                        vt = 1/temp;
                        ae_v_muld(&v.ptr.p_double[1], 1, ae_v_len(1,nv), vt);
                    }
                }
            }
            
            /*
             * Multiple-shift QR step
             */
            for(k=l; k<=i-1; k++)
            {
                
                /*
                 * The first iteration of this loop determines a reflection G
                 * from the vector V and applies it from left and right to H,
                 * thus creating a nonzero bulge below the subdiagonal.
                 *
                 * Each subsequent iteration determines a reflection G to
                 * restore the Hessenberg form in the (K-1)th column, and thus
                 * chases the bulge one step toward the bottom of the active
                 * submatrix. NR is the order of G.
                 */
                nr = ae_minint(ns+1, i-k+1, _state);
                if( k>l )
                {
                    p1 = k-1;
                    p2 = k+nr-1;
                    ae_v_move(&v.ptr.p_double[1], 1, &h->ptr.pp_double[k][p1], h->stride, ae_v_len(1,nr));
                    touchint(&p2, _state);
                }
                generatereflection(&v, nr, &tau, _state);
                if( k>l )
                {
                    h->ptr.pp_double[k][k-1] = v.ptr.p_double[1];
                    for(ii=k+1; ii<=i; ii++)
                    {
                        h->ptr.pp_double[ii][k-1] = (double)(0);
                    }
                }
                v.ptr.p_double[1] = (double)(1);
                
                /*
                 * Apply G from the left to transform the rows of the matrix in
                 * columns K to I2.
                 */
                applyreflectionfromtheleft(h, tau, &v, k, k+nr-1, k, i2, &work, _state);
                
                /*
                 * Apply G from the right to transform the columns of the
                 * matrix in rows I1 to min(K+NR,I).
                 */
                applyreflectionfromtheright(h, tau, &v, i1, ae_minint(k+nr, i, _state), k, k+nr-1, &work, _state);
                if( wantz )
                {
                    
                    /*
                     * Accumulate transformations in the matrix Z
                     */
                    applyreflectionfromtheright(z, tau, &v, 1, n, k, k+nr-1, &work, _state);
                }
            }
        }
        
        /*
         * Failure to converge in remaining number of iterations
         */
        if( failflag )
        {
            *info = i;
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * A submatrix of order <= MAXB in rows and columns L to I has split
         * off. Use the double-shift QR algorithm to handle it.
         */
        hsschur_internalauxschur(wantt, wantz, n, l, i, h, wr, wi, 1, n, z, &work, &workv3, &workc1, &works1, info, _state);
        if( *info>0 )
        {
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * Decrement number of remaining iterations, and return to start of
         * the main loop with a new value of I.
         */
        itn = itn-its;
        i = l-1;
    }
    ae_frame_leave(_state);
}


static void hsschur_internalauxschur(ae_bool wantt,
     ae_bool wantz,
     ae_int_t n,
     ae_int_t ilo,
     ae_int_t ihi,
     /* Real    */ ae_matrix* h,
     /* Real    */ ae_vector* wr,
     /* Real    */ ae_vector* wi,
     ae_int_t iloz,
     ae_int_t ihiz,
     /* Real    */ ae_matrix* z,
     /* Real    */ ae_vector* work,
     /* Real    */ ae_vector* workv3,
     /* Real    */ ae_vector* workc1,
     /* Real    */ ae_vector* works1,
     ae_int_t* info,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t itn;
    ae_int_t its;
    ae_int_t j;
    ae_int_t k;
    ae_int_t l;
    ae_int_t m;
    ae_int_t nh;
    ae_int_t nr;
    ae_int_t nz;
    double ave;
    double cs;
    double disc;
    double h00;
    double h10;
    double h11;
    double h12;
    double h21;
    double h22;
    double h33;
    double h33s;
    double h43h34;
    double h44;
    double h44s;
    double s;
    double smlnum;
    double sn;
    double sum;
    double t1;
    double t2;
    double t3;
    double tst1;
    double unfl;
    double v1;
    double v2;
    double v3;
    ae_bool failflag;
    double dat1;
    double dat2;
    ae_int_t p1;
    double him1im1;
    double him1i;
    double hiim1;
    double hii;
    double wrim1;
    double wri;
    double wiim1;
    double wii;
    double ulp;

    *info = 0;

    *info = 0;
    dat1 = 0.75;
    dat2 = -0.4375;
    ulp = ae_machineepsilon;
    
    /*
     * Quick return if possible
     */
    if( n==0 )
    {
        return;
    }
    if( ilo==ihi )
    {
        wr->ptr.p_double[ilo] = h->ptr.pp_double[ilo][ilo];
        wi->ptr.p_double[ilo] = (double)(0);
        return;
    }
    nh = ihi-ilo+1;
    nz = ihiz-iloz+1;
    
    /*
     * Set machine-dependent constants for the stopping criterion.
     * If norm(H) <= sqrt(MaxRealNumber), overflow should not occur.
     */
    unfl = ae_minrealnumber;
    smlnum = unfl*(nh/ulp);
    
    /*
     * I1 and I2 are the indices of the first row and last column of H
     * to which transformations must be applied. If eigenvalues only are
     * being computed, I1 and I2 are set inside the main loop.
     */
    i1 = 1;
    i2 = n;
    
    /*
     * ITN is the total number of QR iterations allowed.
     */
    itn = 30*nh;
    
    /*
     * The main loop begins here. I is the loop index and decreases from
     * IHI to ILO in steps of 1 or 2. Each iteration of the loop works
     * with the active submatrix in rows and columns L to I.
     * Eigenvalues I+1 to IHI have already converged. Either L = ILO or
     * H(L,L-1) is negligible so that the matrix splits.
     */
    i = ihi;
    for(;;)
    {
        l = ilo;
        if( i<ilo )
        {
            return;
        }
        
        /*
         * Perform QR iterations on rows and columns ILO to I until a
         * submatrix of order 1 or 2 splits off at the bottom because a
         * subdiagonal element has become negligible.
         */
        failflag = ae_true;
        for(its=0; its<=itn; its++)
        {
            
            /*
             * Look for a single small subdiagonal element.
             */
            for(k=i; k>=l+1; k--)
            {
                tst1 = ae_fabs(h->ptr.pp_double[k-1][k-1], _state)+ae_fabs(h->ptr.pp_double[k][k], _state);
                if( ae_fp_eq(tst1,(double)(0)) )
                {
                    tst1 = upperhessenberg1norm(h, l, i, l, i, work, _state);
                }
                if( ae_fp_less_eq(ae_fabs(h->ptr.pp_double[k][k-1], _state),ae_maxreal(ulp*tst1, smlnum, _state)) )
                {
                    break;
                }
            }
            l = k;
            if( l>ilo )
            {
                
                /*
                 * H(L,L-1) is negligible
                 */
                h->ptr.pp_double[l][l-1] = (double)(0);
            }
            
            /*
             * Exit from loop if a submatrix of order 1 or 2 has split off.
             */
            if( l>=i-1 )
            {
                failflag = ae_false;
                break;
            }
            
            /*
             * Now the active submatrix is in rows and columns L to I. If
             * eigenvalues only are being computed, only the active submatrix
             * need be transformed.
             */
            if( its==10||its==20 )
            {
                
                /*
                 * Exceptional shift.
                 */
                s = ae_fabs(h->ptr.pp_double[i][i-1], _state)+ae_fabs(h->ptr.pp_double[i-1][i-2], _state);
                h44 = dat1*s+h->ptr.pp_double[i][i];
                h33 = h44;
                h43h34 = dat2*s*s;
            }
            else
            {
                
                /*
                 * Prepare to use Francis' double shift
                 * (i.e. 2nd degree generalized Rayleigh quotient)
                 */
                h44 = h->ptr.pp_double[i][i];
                h33 = h->ptr.pp_double[i-1][i-1];
                h43h34 = h->ptr.pp_double[i][i-1]*h->ptr.pp_double[i-1][i];
                s = h->ptr.pp_double[i-1][i-2]*h->ptr.pp_double[i-1][i-2];
                disc = (h33-h44)*0.5;
                disc = disc*disc+h43h34;
                if( ae_fp_greater(disc,(double)(0)) )
                {
                    
                    /*
                     * Real roots: use Wilkinson's shift twice
                     */
                    disc = ae_sqrt(disc, _state);
                    ave = 0.5*(h33+h44);
                    if( ae_fp_greater(ae_fabs(h33, _state)-ae_fabs(h44, _state),(double)(0)) )
                    {
                        h33 = h33*h44-h43h34;
                        h44 = h33/(hsschur_extschursign(disc, ave, _state)+ave);
                    }
                    else
                    {
                        h44 = hsschur_extschursign(disc, ave, _state)+ave;
                    }
                    h33 = h44;
                    h43h34 = (double)(0);
                }
            }
            
            /*
             * Look for two consecutive small subdiagonal elements.
             */
            for(m=i-2; m>=l; m--)
            {
                
                /*
                 * Determine the effect of starting the double-shift QR
                 * iteration at row M, and see if this would make H(M,M-1)
                 * negligible.
                 */
                h11 = h->ptr.pp_double[m][m];
                h22 = h->ptr.pp_double[m+1][m+1];
                h21 = h->ptr.pp_double[m+1][m];
                h12 = h->ptr.pp_double[m][m+1];
                h44s = h44-h11;
                h33s = h33-h11;
                v1 = (h33s*h44s-h43h34)/h21+h12;
                v2 = h22-h11-h33s-h44s;
                v3 = h->ptr.pp_double[m+2][m+1];
                s = ae_fabs(v1, _state)+ae_fabs(v2, _state)+ae_fabs(v3, _state);
                v1 = v1/s;
                v2 = v2/s;
                v3 = v3/s;
                workv3->ptr.p_double[1] = v1;
                workv3->ptr.p_double[2] = v2;
                workv3->ptr.p_double[3] = v3;
                if( m==l )
                {
                    break;
                }
                h00 = h->ptr.pp_double[m-1][m-1];
                h10 = h->ptr.pp_double[m][m-1];
                tst1 = ae_fabs(v1, _state)*(ae_fabs(h00, _state)+ae_fabs(h11, _state)+ae_fabs(h22, _state));
                if( ae_fp_less_eq(ae_fabs(h10, _state)*(ae_fabs(v2, _state)+ae_fabs(v3, _state)),ulp*tst1) )
                {
                    break;
                }
            }
            
            /*
             * Double-shift QR step
             */
            for(k=m; k<=i-1; k++)
            {
                
                /*
                 * The first iteration of this loop determines a reflection G
                 * from the vector V and applies it from left and right to H,
                 * thus creating a nonzero bulge below the subdiagonal.
                 *
                 * Each subsequent iteration determines a reflection G to
                 * restore the Hessenberg form in the (K-1)th column, and thus
                 * chases the bulge one step toward the bottom of the active
                 * submatrix. NR is the order of G.
                 */
                nr = ae_minint(3, i-k+1, _state);
                if( k>m )
                {
                    for(p1=1; p1<=nr; p1++)
                    {
                        workv3->ptr.p_double[p1] = h->ptr.pp_double[k+p1-1][k-1];
                    }
                }
                generatereflection(workv3, nr, &t1, _state);
                if( k>m )
                {
                    h->ptr.pp_double[k][k-1] = workv3->ptr.p_double[1];
                    h->ptr.pp_double[k+1][k-1] = (double)(0);
                    if( k<i-1 )
                    {
                        h->ptr.pp_double[k+2][k-1] = (double)(0);
                    }
                }
                else
                {
                    if( m>l )
                    {
                        h->ptr.pp_double[k][k-1] = -h->ptr.pp_double[k][k-1];
                    }
                }
                v2 = workv3->ptr.p_double[2];
                t2 = t1*v2;
                if( nr==3 )
                {
                    v3 = workv3->ptr.p_double[3];
                    t3 = t1*v3;
                    
                    /*
                     * Apply G from the left to transform the rows of the matrix
                     * in columns K to I2.
                     */
                    for(j=k; j<=i2; j++)
                    {
                        sum = h->ptr.pp_double[k][j]+v2*h->ptr.pp_double[k+1][j]+v3*h->ptr.pp_double[k+2][j];
                        h->ptr.pp_double[k][j] = h->ptr.pp_double[k][j]-sum*t1;
                        h->ptr.pp_double[k+1][j] = h->ptr.pp_double[k+1][j]-sum*t2;
                        h->ptr.pp_double[k+2][j] = h->ptr.pp_double[k+2][j]-sum*t3;
                    }
                    
                    /*
                     * Apply G from the right to transform the columns of the
                     * matrix in rows I1 to min(K+3,I).
                     */
                    for(j=i1; j<=ae_minint(k+3, i, _state); j++)
                    {
                        sum = h->ptr.pp_double[j][k]+v2*h->ptr.pp_double[j][k+1]+v3*h->ptr.pp_double[j][k+2];
                        h->ptr.pp_double[j][k] = h->ptr.pp_double[j][k]-sum*t1;
                        h->ptr.pp_double[j][k+1] = h->ptr.pp_double[j][k+1]-sum*t2;
                        h->ptr.pp_double[j][k+2] = h->ptr.pp_double[j][k+2]-sum*t3;
                    }
                    if( wantz )
                    {
                        
                        /*
                         * Accumulate transformations in the matrix Z
                         */
                        for(j=iloz; j<=ihiz; j++)
                        {
                            sum = z->ptr.pp_double[j][k]+v2*z->ptr.pp_double[j][k+1]+v3*z->ptr.pp_double[j][k+2];
                            z->ptr.pp_double[j][k] = z->ptr.pp_double[j][k]-sum*t1;
                            z->ptr.pp_double[j][k+1] = z->ptr.pp_double[j][k+1]-sum*t2;
                            z->ptr.pp_double[j][k+2] = z->ptr.pp_double[j][k+2]-sum*t3;
                        }
                    }
                }
                else
                {
                    if( nr==2 )
                    {
                        
                        /*
                         * Apply G from the left to transform the rows of the matrix
                         * in columns K to I2.
                         */
                        for(j=k; j<=i2; j++)
                        {
                            sum = h->ptr.pp_double[k][j]+v2*h->ptr.pp_double[k+1][j];
                            h->ptr.pp_double[k][j] = h->ptr.pp_double[k][j]-sum*t1;
                            h->ptr.pp_double[k+1][j] = h->ptr.pp_double[k+1][j]-sum*t2;
                        }
                        
                        /*
                         * Apply G from the right to transform the columns of the
                         * matrix in rows I1 to min(K+3,I).
                         */
                        for(j=i1; j<=i; j++)
                        {
                            sum = h->ptr.pp_double[j][k]+v2*h->ptr.pp_double[j][k+1];
                            h->ptr.pp_double[j][k] = h->ptr.pp_double[j][k]-sum*t1;
                            h->ptr.pp_double[j][k+1] = h->ptr.pp_double[j][k+1]-sum*t2;
                        }
                        if( wantz )
                        {
                            
                            /*
                             * Accumulate transformations in the matrix Z
                             */
                            for(j=iloz; j<=ihiz; j++)
                            {
                                sum = z->ptr.pp_double[j][k]+v2*z->ptr.pp_double[j][k+1];
                                z->ptr.pp_double[j][k] = z->ptr.pp_double[j][k]-sum*t1;
                                z->ptr.pp_double[j][k+1] = z->ptr.pp_double[j][k+1]-sum*t2;
                            }
                        }
                    }
                }
            }
        }
        if( failflag )
        {
            
            /*
             * Failure to converge in remaining number of iterations
             */
            *info = i;
            return;
        }
        if( l==i )
        {
            
            /*
             * H(I,I-1) is negligible: one eigenvalue has converged.
             */
            wr->ptr.p_double[i] = h->ptr.pp_double[i][i];
            wi->ptr.p_double[i] = (double)(0);
        }
        else
        {
            if( l==i-1 )
            {
                
                /*
                 * H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
                 *
                 *        Transform the 2-by-2 submatrix to standard Schur form,
                 *        and compute and store the eigenvalues.
                 */
                him1im1 = h->ptr.pp_double[i-1][i-1];
                him1i = h->ptr.pp_double[i-1][i];
                hiim1 = h->ptr.pp_double[i][i-1];
                hii = h->ptr.pp_double[i][i];
                hsschur_aux2x2schur(&him1im1, &him1i, &hiim1, &hii, &wrim1, &wiim1, &wri, &wii, &cs, &sn, _state);
                wr->ptr.p_double[i-1] = wrim1;
                wi->ptr.p_double[i-1] = wiim1;
                wr->ptr.p_double[i] = wri;
                wi->ptr.p_double[i] = wii;
                h->ptr.pp_double[i-1][i-1] = him1im1;
                h->ptr.pp_double[i-1][i] = him1i;
                h->ptr.pp_double[i][i-1] = hiim1;
                h->ptr.pp_double[i][i] = hii;
                if( wantt )
                {
                    
                    /*
                     * Apply the transformation to the rest of H.
                     */
                    if( i2>i )
                    {
                        workc1->ptr.p_double[1] = cs;
                        works1->ptr.p_double[1] = sn;
                        applyrotationsfromtheleft(ae_true, i-1, i, i+1, i2, workc1, works1, h, work, _state);
                    }
                    workc1->ptr.p_double[1] = cs;
                    works1->ptr.p_double[1] = sn;
                    applyrotationsfromtheright(ae_true, i1, i-2, i-1, i, workc1, works1, h, work, _state);
                }
                if( wantz )
                {
                    
                    /*
                     * Apply the transformation to Z.
                     */
                    workc1->ptr.p_double[1] = cs;
                    works1->ptr.p_double[1] = sn;
                    applyrotationsfromtheright(ae_true, iloz, iloz+nz-1, i-1, i, workc1, works1, z, work, _state);
                }
            }
        }
        
        /*
         * Decrement number of remaining iterations, and return to start of
         * the main loop with new value of I.
         */
        itn = itn-its;
        i = l-1;
    }
}


static void hsschur_aux2x2schur(double* a,
     double* b,
     double* c,
     double* d,
     double* rt1r,
     double* rt1i,
     double* rt2r,
     double* rt2i,
     double* cs,
     double* sn,
     ae_state *_state)
{
    double multpl;
    double aa;
    double bb;
    double bcmax;
    double bcmis;
    double cc;
    double cs1;
    double dd;
    double eps;
    double p;
    double sab;
    double sac;
    double scl;
    double sigma;
    double sn1;
    double tau;
    double temp;
    double z;

    *rt1r = 0;
    *rt1i = 0;
    *rt2r = 0;
    *rt2i = 0;
    *cs = 0;
    *sn = 0;

    multpl = 4.0;
    eps = ae_machineepsilon;
    if( ae_fp_eq(*c,(double)(0)) )
    {
        *cs = (double)(1);
        *sn = (double)(0);
    }
    else
    {
        if( ae_fp_eq(*b,(double)(0)) )
        {
            
            /*
             * Swap rows and columns
             */
            *cs = (double)(0);
            *sn = (double)(1);
            temp = *d;
            *d = *a;
            *a = temp;
            *b = -*c;
            *c = (double)(0);
        }
        else
        {
            if( ae_fp_eq(*a-(*d),(double)(0))&&hsschur_extschursigntoone(*b, _state)!=hsschur_extschursigntoone(*c, _state) )
            {
                *cs = (double)(1);
                *sn = (double)(0);
            }
            else
            {
                temp = *a-(*d);
                p = 0.5*temp;
                bcmax = ae_maxreal(ae_fabs(*b, _state), ae_fabs(*c, _state), _state);
                bcmis = ae_minreal(ae_fabs(*b, _state), ae_fabs(*c, _state), _state)*hsschur_extschursigntoone(*b, _state)*hsschur_extschursigntoone(*c, _state);
                scl = ae_maxreal(ae_fabs(p, _state), bcmax, _state);
                z = p/scl*p+bcmax/scl*bcmis;
                
                /*
                 * If Z is of the order of the machine accuracy, postpone the
                 * decision on the nature of eigenvalues
                 */
                if( ae_fp_greater_eq(z,multpl*eps) )
                {
                    
                    /*
                     * Real eigenvalues. Compute A and D.
                     */
                    z = p+hsschur_extschursign(ae_sqrt(scl, _state)*ae_sqrt(z, _state), p, _state);
                    *a = *d+z;
                    *d = *d-bcmax/z*bcmis;
                    
                    /*
                     * Compute B and the rotation matrix
                     */
                    tau = pythag2(*c, z, _state);
                    *cs = z/tau;
                    *sn = *c/tau;
                    *b = *b-(*c);
                    *c = (double)(0);
                }
                else
                {
                    
                    /*
                     * Complex eigenvalues, or real (almost) equal eigenvalues.
                     * Make diagonal elements equal.
                     */
                    sigma = *b+(*c);
                    tau = pythag2(sigma, temp, _state);
                    *cs = ae_sqrt(0.5*(1+ae_fabs(sigma, _state)/tau), _state);
                    *sn = -p/(tau*(*cs))*hsschur_extschursign((double)(1), sigma, _state);
                    
                    /*
                     * Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
                     *         [ CC  DD ]   [ C  D ] [ SN  CS ]
                     */
                    aa = *a*(*cs)+*b*(*sn);
                    bb = -*a*(*sn)+*b*(*cs);
                    cc = *c*(*cs)+*d*(*sn);
                    dd = -*c*(*sn)+*d*(*cs);
                    
                    /*
                     * Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
                     *         [ C  D ]   [-SN  CS ] [ CC  DD ]
                     */
                    *a = aa*(*cs)+cc*(*sn);
                    *b = bb*(*cs)+dd*(*sn);
                    *c = -aa*(*sn)+cc*(*cs);
                    *d = -bb*(*sn)+dd*(*cs);
                    temp = 0.5*(*a+(*d));
                    *a = temp;
                    *d = temp;
                    if( ae_fp_neq(*c,(double)(0)) )
                    {
                        if( ae_fp_neq(*b,(double)(0)) )
                        {
                            if( hsschur_extschursigntoone(*b, _state)==hsschur_extschursigntoone(*c, _state) )
                            {
                                
                                /*
                                 * Real eigenvalues: reduce to upper triangular form
                                 */
                                sab = ae_sqrt(ae_fabs(*b, _state), _state);
                                sac = ae_sqrt(ae_fabs(*c, _state), _state);
                                p = hsschur_extschursign(sab*sac, *c, _state);
                                tau = 1/ae_sqrt(ae_fabs(*b+(*c), _state), _state);
                                *a = temp+p;
                                *d = temp-p;
                                *b = *b-(*c);
                                *c = (double)(0);
                                cs1 = sab*tau;
                                sn1 = sac*tau;
                                temp = *cs*cs1-*sn*sn1;
                                *sn = *cs*sn1+*sn*cs1;
                                *cs = temp;
                            }
                        }
                        else
                        {
                            *b = -*c;
                            *c = (double)(0);
                            temp = *cs;
                            *cs = -*sn;
                            *sn = temp;
                        }
                    }
                }
            }
        }
    }
    
    /*
     * Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
     */
    *rt1r = *a;
    *rt2r = *d;
    if( ae_fp_eq(*c,(double)(0)) )
    {
        *rt1i = (double)(0);
        *rt2i = (double)(0);
    }
    else
    {
        *rt1i = ae_sqrt(ae_fabs(*b, _state), _state)*ae_sqrt(ae_fabs(*c, _state), _state);
        *rt2i = -*rt1i;
    }
}


static double hsschur_extschursign(double a, double b, ae_state *_state)
{
    double result;


    if( ae_fp_greater_eq(b,(double)(0)) )
    {
        result = ae_fabs(a, _state);
    }
    else
    {
        result = -ae_fabs(a, _state);
    }
    return result;
}


static ae_int_t hsschur_extschursigntoone(double b, ae_state *_state)
{
    ae_int_t result;


    if( ae_fp_greater_eq(b,(double)(0)) )
    {
        result = 1;
    }
    else
    {
        result = -1;
    }
    return result;
}




/*************************************************************************
Utility subroutine performing the "safe" solution of system of linear
equations with triangular coefficient matrices.

The subroutine uses scaling and solves the scaled system A*x=s*b (where  s
is  a  scalar  value)  instead  of  A*x=b,  choosing  s  so  that x can be
represented by a floating-point number. The closer the system  gets  to  a
singular, the less s is. If the system is singular, s=0 and x contains the
non-trivial solution of equation A*x=0.

The feature of an algorithm is that it could not cause an  overflow  or  a
division by zero regardless of the matrix used as the input.

The algorithm can solve systems of equations with  upper/lower  triangular
matrices,  with/without unit diagonal, and systems of type A*x=b or A'*x=b
(where A' is a transposed matrix A).

Input parameters:
    A       -   system matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    X       -   right-hand member of a system.
                Array whose index ranges within [0..N-1].
    IsUpper -   matrix type. If it is True, the system matrix is the upper
                triangular and is located in  the  corresponding  part  of
                matrix A.
    Trans   -   problem type. If it is True, the problem to be  solved  is
                A'*x=b, otherwise it is A*x=b.
    Isunit  -   matrix type. If it is True, the system matrix has  a  unit
                diagonal (the elements on the main diagonal are  not  used
                in the calculation process), otherwise the matrix is considered
                to be a general triangular matrix.

Output parameters:
    X       -   solution. Array whose index ranges within [0..N-1].
    S       -   scaling factor.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************/
void rmatrixtrsafesolve(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     double* s,
     ae_bool isupper,
     ae_bool istrans,
     ae_bool isunit,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool normin;
    ae_vector cnorm;
    ae_matrix a1;
    ae_vector x1;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    *s = 0;
    ae_vector_init(&cnorm, 0, DT_REAL, _state);
    ae_matrix_init(&a1, 0, 0, DT_REAL, _state);
    ae_vector_init(&x1, 0, DT_REAL, _state);

    
    /*
     * From 0-based to 1-based
     */
    normin = ae_false;
    ae_matrix_set_length(&a1, n+1, n+1, _state);
    ae_vector_set_length(&x1, n+1, _state);
    for(i=1; i<=n; i++)
    {
        ae_v_move(&a1.ptr.pp_double[i][1], 1, &a->ptr.pp_double[i-1][0], 1, ae_v_len(1,n));
    }
    ae_v_move(&x1.ptr.p_double[1], 1, &x->ptr.p_double[0], 1, ae_v_len(1,n));
    
    /*
     * Solve 1-based
     */
    safesolvetriangular(&a1, n, &x1, s, isupper, istrans, isunit, normin, &cnorm, _state);
    
    /*
     * From 1-based to 0-based
     */
    ae_v_move(&x->ptr.p_double[0], 1, &x1.ptr.p_double[1], 1, ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Obsolete 1-based subroutine.
See RMatrixTRSafeSolve for 0-based replacement.
*************************************************************************/
void safesolvetriangular(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     double* s,
     ae_bool isupper,
     ae_bool istrans,
     ae_bool isunit,
     ae_bool normin,
     /* Real    */ ae_vector* cnorm,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t imax;
    ae_int_t j;
    ae_int_t jfirst;
    ae_int_t jinc;
    ae_int_t jlast;
    ae_int_t jm1;
    ae_int_t jp1;
    ae_int_t ip1;
    ae_int_t im1;
    ae_int_t k;
    ae_int_t flg;
    double v;
    double vd;
    double bignum;
    double grow;
    double rec;
    double smlnum;
    double sumj;
    double tjj;
    double tjjs;
    double tmax;
    double tscal;
    double uscal;
    double xbnd;
    double xj;
    double xmax;
    ae_bool notran;
    ae_bool upper;
    ae_bool nounit;

    *s = 0;

    upper = isupper;
    notran = !istrans;
    nounit = !isunit;
    
    /*
     * these initializers are not really necessary,
     * but without them compiler complains about uninitialized locals
     */
    tjjs = (double)(0);
    
    /*
     * Quick return if possible
     */
    if( n==0 )
    {
        return;
    }
    
    /*
     * Determine machine dependent parameters to control overflow.
     */
    smlnum = ae_minrealnumber/(ae_machineepsilon*2);
    bignum = 1/smlnum;
    *s = (double)(1);
    if( !normin )
    {
        ae_vector_set_length(cnorm, n+1, _state);
        
        /*
         * Compute the 1-norm of each column, not including the diagonal.
         */
        if( upper )
        {
            
            /*
             * A is upper triangular.
             */
            for(j=1; j<=n; j++)
            {
                v = (double)(0);
                for(k=1; k<=j-1; k++)
                {
                    v = v+ae_fabs(a->ptr.pp_double[k][j], _state);
                }
                cnorm->ptr.p_double[j] = v;
            }
        }
        else
        {
            
            /*
             * A is lower triangular.
             */
            for(j=1; j<=n-1; j++)
            {
                v = (double)(0);
                for(k=j+1; k<=n; k++)
                {
                    v = v+ae_fabs(a->ptr.pp_double[k][j], _state);
                }
                cnorm->ptr.p_double[j] = v;
            }
            cnorm->ptr.p_double[n] = (double)(0);
        }
    }
    
    /*
     * Scale the column norms by TSCAL if the maximum element in CNORM is
     * greater than BIGNUM.
     */
    imax = 1;
    for(k=2; k<=n; k++)
    {
        if( ae_fp_greater(cnorm->ptr.p_double[k],cnorm->ptr.p_double[imax]) )
        {
            imax = k;
        }
    }
    tmax = cnorm->ptr.p_double[imax];
    if( ae_fp_less_eq(tmax,bignum) )
    {
        tscal = (double)(1);
    }
    else
    {
        tscal = 1/(smlnum*tmax);
        ae_v_muld(&cnorm->ptr.p_double[1], 1, ae_v_len(1,n), tscal);
    }
    
    /*
     * Compute a bound on the computed solution vector to see if the
     * Level 2 BLAS routine DTRSV can be used.
     */
    j = 1;
    for(k=2; k<=n; k++)
    {
        if( ae_fp_greater(ae_fabs(x->ptr.p_double[k], _state),ae_fabs(x->ptr.p_double[j], _state)) )
        {
            j = k;
        }
    }
    xmax = ae_fabs(x->ptr.p_double[j], _state);
    xbnd = xmax;
    if( notran )
    {
        
        /*
         * Compute the growth in A * x = b.
         */
        if( upper )
        {
            jfirst = n;
            jlast = 1;
            jinc = -1;
        }
        else
        {
            jfirst = 1;
            jlast = n;
            jinc = 1;
        }
        if( ae_fp_neq(tscal,(double)(1)) )
        {
            grow = (double)(0);
        }
        else
        {
            if( nounit )
            {
                
                /*
                 * A is non-unit triangular.
                 *
                 * Compute GROW = 1/G(j) and XBND = 1/M(j).
                 * Initially, G(0) = max{x(i), i=1,...,n}.
                 */
                grow = 1/ae_maxreal(xbnd, smlnum, _state);
                xbnd = grow;
                j = jfirst;
                while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
                {
                    
                    /*
                     * Exit the loop if the growth factor is too small.
                     */
                    if( ae_fp_less_eq(grow,smlnum) )
                    {
                        break;
                    }
                    
                    /*
                     * M(j) = G(j-1) / abs(A(j,j))
                     */
                    tjj = ae_fabs(a->ptr.pp_double[j][j], _state);
                    xbnd = ae_minreal(xbnd, ae_minreal((double)(1), tjj, _state)*grow, _state);
                    if( ae_fp_greater_eq(tjj+cnorm->ptr.p_double[j],smlnum) )
                    {
                        
                        /*
                         * G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
                         */
                        grow = grow*(tjj/(tjj+cnorm->ptr.p_double[j]));
                    }
                    else
                    {
                        
                        /*
                         * G(j) could overflow, set GROW to 0.
                         */
                        grow = (double)(0);
                    }
                    if( j==jlast )
                    {
                        grow = xbnd;
                    }
                    j = j+jinc;
                }
            }
            else
            {
                
                /*
                 * A is unit triangular.
                 *
                 * Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
                 */
                grow = ae_minreal((double)(1), 1/ae_maxreal(xbnd, smlnum, _state), _state);
                j = jfirst;
                while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
                {
                    
                    /*
                     * Exit the loop if the growth factor is too small.
                     */
                    if( ae_fp_less_eq(grow,smlnum) )
                    {
                        break;
                    }
                    
                    /*
                     * G(j) = G(j-1)*( 1 + CNORM(j) )
                     */
                    grow = grow*(1/(1+cnorm->ptr.p_double[j]));
                    j = j+jinc;
                }
            }
        }
    }
    else
    {
        
        /*
         * Compute the growth in A' * x = b.
         */
        if( upper )
        {
            jfirst = 1;
            jlast = n;
            jinc = 1;
        }
        else
        {
            jfirst = n;
            jlast = 1;
            jinc = -1;
        }
        if( ae_fp_neq(tscal,(double)(1)) )
        {
            grow = (double)(0);
        }
        else
        {
            if( nounit )
            {
                
                /*
                 * A is non-unit triangular.
                 *
                 * Compute GROW = 1/G(j) and XBND = 1/M(j).
                 * Initially, M(0) = max{x(i), i=1,...,n}.
                 */
                grow = 1/ae_maxreal(xbnd, smlnum, _state);
                xbnd = grow;
                j = jfirst;
                while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
                {
                    
                    /*
                     * Exit the loop if the growth factor is too small.
                     */
                    if( ae_fp_less_eq(grow,smlnum) )
                    {
                        break;
                    }
                    
                    /*
                     * G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
                     */
                    xj = 1+cnorm->ptr.p_double[j];
                    grow = ae_minreal(grow, xbnd/xj, _state);
                    
                    /*
                     * M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
                     */
                    tjj = ae_fabs(a->ptr.pp_double[j][j], _state);
                    if( ae_fp_greater(xj,tjj) )
                    {
                        xbnd = xbnd*(tjj/xj);
                    }
                    if( j==jlast )
                    {
                        grow = ae_minreal(grow, xbnd, _state);
                    }
                    j = j+jinc;
                }
            }
            else
            {
                
                /*
                 * A is unit triangular.
                 *
                 * Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
                 */
                grow = ae_minreal((double)(1), 1/ae_maxreal(xbnd, smlnum, _state), _state);
                j = jfirst;
                while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
                {
                    
                    /*
                     * Exit the loop if the growth factor is too small.
                     */
                    if( ae_fp_less_eq(grow,smlnum) )
                    {
                        break;
                    }
                    
                    /*
                     * G(j) = ( 1 + CNORM(j) )*G(j-1)
                     */
                    xj = 1+cnorm->ptr.p_double[j];
                    grow = grow/xj;
                    j = j+jinc;
                }
            }
        }
    }
    if( ae_fp_greater(grow*tscal,smlnum) )
    {
        
        /*
         * Use the Level 2 BLAS solve if the reciprocal of the bound on
         * elements of X is not too small.
         */
        if( (upper&&notran)||(!upper&&!notran) )
        {
            if( nounit )
            {
                vd = a->ptr.pp_double[n][n];
            }
            else
            {
                vd = (double)(1);
            }
            x->ptr.p_double[n] = x->ptr.p_double[n]/vd;
            for(i=n-1; i>=1; i--)
            {
                ip1 = i+1;
                if( upper )
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[i][ip1], 1, &x->ptr.p_double[ip1], 1, ae_v_len(ip1,n));
                }
                else
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[ip1][i], a->stride, &x->ptr.p_double[ip1], 1, ae_v_len(ip1,n));
                }
                if( nounit )
                {
                    vd = a->ptr.pp_double[i][i];
                }
                else
                {
                    vd = (double)(1);
                }
                x->ptr.p_double[i] = (x->ptr.p_double[i]-v)/vd;
            }
        }
        else
        {
            if( nounit )
            {
                vd = a->ptr.pp_double[1][1];
            }
            else
            {
                vd = (double)(1);
            }
            x->ptr.p_double[1] = x->ptr.p_double[1]/vd;
            for(i=2; i<=n; i++)
            {
                im1 = i-1;
                if( upper )
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[1][i], a->stride, &x->ptr.p_double[1], 1, ae_v_len(1,im1));
                }
                else
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[i][1], 1, &x->ptr.p_double[1], 1, ae_v_len(1,im1));
                }
                if( nounit )
                {
                    vd = a->ptr.pp_double[i][i];
                }
                else
                {
                    vd = (double)(1);
                }
                x->ptr.p_double[i] = (x->ptr.p_double[i]-v)/vd;
            }
        }
    }
    else
    {
        
        /*
         * Use a Level 1 BLAS solve, scaling intermediate results.
         */
        if( ae_fp_greater(xmax,bignum) )
        {
            
            /*
             * Scale X so that its components are less than or equal to
             * BIGNUM in absolute value.
             */
            *s = bignum/xmax;
            ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), *s);
            xmax = bignum;
        }
        if( notran )
        {
            
            /*
             * Solve A * x = b
             */
            j = jfirst;
            while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
            {
                
                /*
                 * Compute x(j) = b(j) / A(j,j), scaling x if necessary.
                 */
                xj = ae_fabs(x->ptr.p_double[j], _state);
                flg = 0;
                if( nounit )
                {
                    tjjs = a->ptr.pp_double[j][j]*tscal;
                }
                else
                {
                    tjjs = tscal;
                    if( ae_fp_eq(tscal,(double)(1)) )
                    {
                        flg = 100;
                    }
                }
                if( flg!=100 )
                {
                    tjj = ae_fabs(tjjs, _state);
                    if( ae_fp_greater(tjj,smlnum) )
                    {
                        
                        /*
                         * abs(A(j,j)) > SMLNUM:
                         */
                        if( ae_fp_less(tjj,(double)(1)) )
                        {
                            if( ae_fp_greater(xj,tjj*bignum) )
                            {
                                
                                /*
                                 * Scale x by 1/b(j).
                                 */
                                rec = 1/xj;
                                ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                                *s = *s*rec;
                                xmax = xmax*rec;
                            }
                        }
                        x->ptr.p_double[j] = x->ptr.p_double[j]/tjjs;
                        xj = ae_fabs(x->ptr.p_double[j], _state);
                    }
                    else
                    {
                        if( ae_fp_greater(tjj,(double)(0)) )
                        {
                            
                            /*
                             * 0 < abs(A(j,j)) <= SMLNUM:
                             */
                            if( ae_fp_greater(xj,tjj*bignum) )
                            {
                                
                                /*
                                 * Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                                 * to avoid overflow when dividing by A(j,j).
                                 */
                                rec = tjj*bignum/xj;
                                if( ae_fp_greater(cnorm->ptr.p_double[j],(double)(1)) )
                                {
                                    
                                    /*
                                     * Scale by 1/CNORM(j) to avoid overflow when
                                     * multiplying x(j) times column j.
                                     */
                                    rec = rec/cnorm->ptr.p_double[j];
                                }
                                ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                                *s = *s*rec;
                                xmax = xmax*rec;
                            }
                            x->ptr.p_double[j] = x->ptr.p_double[j]/tjjs;
                            xj = ae_fabs(x->ptr.p_double[j], _state);
                        }
                        else
                        {
                            
                            /*
                             * A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                             * scale = 0, and compute a solution to A*x = 0.
                             */
                            for(i=1; i<=n; i++)
                            {
                                x->ptr.p_double[i] = (double)(0);
                            }
                            x->ptr.p_double[j] = (double)(1);
                            xj = (double)(1);
                            *s = (double)(0);
                            xmax = (double)(0);
                        }
                    }
                }
                
                /*
                 * Scale x if necessary to avoid overflow when adding a
                 * multiple of column j of A.
                 */
                if( ae_fp_greater(xj,(double)(1)) )
                {
                    rec = 1/xj;
                    if( ae_fp_greater(cnorm->ptr.p_double[j],(bignum-xmax)*rec) )
                    {
                        
                        /*
                         * Scale x by 1/(2*abs(x(j))).
                         */
                        rec = rec*0.5;
                        ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                        *s = *s*rec;
                    }
                }
                else
                {
                    if( ae_fp_greater(xj*cnorm->ptr.p_double[j],bignum-xmax) )
                    {
                        
                        /*
                         * Scale x by 1/2.
                         */
                        ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), 0.5);
                        *s = *s*0.5;
                    }
                }
                if( upper )
                {
                    if( j>1 )
                    {
                        
                        /*
                         * Compute the update
                         * x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
                         */
                        v = x->ptr.p_double[j]*tscal;
                        jm1 = j-1;
                        ae_v_subd(&x->ptr.p_double[1], 1, &a->ptr.pp_double[1][j], a->stride, ae_v_len(1,jm1), v);
                        i = 1;
                        for(k=2; k<=j-1; k++)
                        {
                            if( ae_fp_greater(ae_fabs(x->ptr.p_double[k], _state),ae_fabs(x->ptr.p_double[i], _state)) )
                            {
                                i = k;
                            }
                        }
                        xmax = ae_fabs(x->ptr.p_double[i], _state);
                    }
                }
                else
                {
                    if( j<n )
                    {
                        
                        /*
                         * Compute the update
                         * x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
                         */
                        jp1 = j+1;
                        v = x->ptr.p_double[j]*tscal;
                        ae_v_subd(&x->ptr.p_double[jp1], 1, &a->ptr.pp_double[jp1][j], a->stride, ae_v_len(jp1,n), v);
                        i = j+1;
                        for(k=j+2; k<=n; k++)
                        {
                            if( ae_fp_greater(ae_fabs(x->ptr.p_double[k], _state),ae_fabs(x->ptr.p_double[i], _state)) )
                            {
                                i = k;
                            }
                        }
                        xmax = ae_fabs(x->ptr.p_double[i], _state);
                    }
                }
                j = j+jinc;
            }
        }
        else
        {
            
            /*
             * Solve A' * x = b
             */
            j = jfirst;
            while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
            {
                
                /*
                 * Compute x(j) = b(j) - sum A(k,j)*x(k).
                 *   k<>j
                 */
                xj = ae_fabs(x->ptr.p_double[j], _state);
                uscal = tscal;
                rec = 1/ae_maxreal(xmax, (double)(1), _state);
                if( ae_fp_greater(cnorm->ptr.p_double[j],(bignum-xj)*rec) )
                {
                    
                    /*
                     * If x(j) could overflow, scale x by 1/(2*XMAX).
                     */
                    rec = rec*0.5;
                    if( nounit )
                    {
                        tjjs = a->ptr.pp_double[j][j]*tscal;
                    }
                    else
                    {
                        tjjs = tscal;
                    }
                    tjj = ae_fabs(tjjs, _state);
                    if( ae_fp_greater(tjj,(double)(1)) )
                    {
                        
                        /*
                         * Divide by A(j,j) when scaling x if A(j,j) > 1.
                         */
                        rec = ae_minreal((double)(1), rec*tjj, _state);
                        uscal = uscal/tjjs;
                    }
                    if( ae_fp_less(rec,(double)(1)) )
                    {
                        ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                        *s = *s*rec;
                        xmax = xmax*rec;
                    }
                }
                sumj = (double)(0);
                if( ae_fp_eq(uscal,(double)(1)) )
                {
                    
                    /*
                     * If the scaling needed for A in the dot product is 1,
                     * call DDOT to perform the dot product.
                     */
                    if( upper )
                    {
                        if( j>1 )
                        {
                            jm1 = j-1;
                            sumj = ae_v_dotproduct(&a->ptr.pp_double[1][j], a->stride, &x->ptr.p_double[1], 1, ae_v_len(1,jm1));
                        }
                        else
                        {
                            sumj = (double)(0);
                        }
                    }
                    else
                    {
                        if( j<n )
                        {
                            jp1 = j+1;
                            sumj = ae_v_dotproduct(&a->ptr.pp_double[jp1][j], a->stride, &x->ptr.p_double[jp1], 1, ae_v_len(jp1,n));
                        }
                    }
                }
                else
                {
                    
                    /*
                     * Otherwise, use in-line code for the dot product.
                     */
                    if( upper )
                    {
                        for(i=1; i<=j-1; i++)
                        {
                            v = a->ptr.pp_double[i][j]*uscal;
                            sumj = sumj+v*x->ptr.p_double[i];
                        }
                    }
                    else
                    {
                        if( j<n )
                        {
                            for(i=j+1; i<=n; i++)
                            {
                                v = a->ptr.pp_double[i][j]*uscal;
                                sumj = sumj+v*x->ptr.p_double[i];
                            }
                        }
                    }
                }
                if( ae_fp_eq(uscal,tscal) )
                {
                    
                    /*
                     * Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
                     * was not used to scale the dotproduct.
                     */
                    x->ptr.p_double[j] = x->ptr.p_double[j]-sumj;
                    xj = ae_fabs(x->ptr.p_double[j], _state);
                    flg = 0;
                    if( nounit )
                    {
                        tjjs = a->ptr.pp_double[j][j]*tscal;
                    }
                    else
                    {
                        tjjs = tscal;
                        if( ae_fp_eq(tscal,(double)(1)) )
                        {
                            flg = 150;
                        }
                    }
                    
                    /*
                     * Compute x(j) = x(j) / A(j,j), scaling if necessary.
                     */
                    if( flg!=150 )
                    {
                        tjj = ae_fabs(tjjs, _state);
                        if( ae_fp_greater(tjj,smlnum) )
                        {
                            
                            /*
                             * abs(A(j,j)) > SMLNUM:
                             */
                            if( ae_fp_less(tjj,(double)(1)) )
                            {
                                if( ae_fp_greater(xj,tjj*bignum) )
                                {
                                    
                                    /*
                                     * Scale X by 1/abs(x(j)).
                                     */
                                    rec = 1/xj;
                                    ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                                    *s = *s*rec;
                                    xmax = xmax*rec;
                                }
                            }
                            x->ptr.p_double[j] = x->ptr.p_double[j]/tjjs;
                        }
                        else
                        {
                            if( ae_fp_greater(tjj,(double)(0)) )
                            {
                                
                                /*
                                 * 0 < abs(A(j,j)) <= SMLNUM:
                                 */
                                if( ae_fp_greater(xj,tjj*bignum) )
                                {
                                    
                                    /*
                                     * Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
                                     */
                                    rec = tjj*bignum/xj;
                                    ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                                    *s = *s*rec;
                                    xmax = xmax*rec;
                                }
                                x->ptr.p_double[j] = x->ptr.p_double[j]/tjjs;
                            }
                            else
                            {
                                
                                /*
                                 * A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                                 * scale = 0, and compute a solution to A'*x = 0.
                                 */
                                for(i=1; i<=n; i++)
                                {
                                    x->ptr.p_double[i] = (double)(0);
                                }
                                x->ptr.p_double[j] = (double)(1);
                                *s = (double)(0);
                                xmax = (double)(0);
                            }
                        }
                    }
                }
                else
                {
                    
                    /*
                     * Compute x(j) := x(j) / A(j,j)  - sumj if the dot
                     * product has already been divided by 1/A(j,j).
                     */
                    x->ptr.p_double[j] = x->ptr.p_double[j]/tjjs-sumj;
                }
                xmax = ae_maxreal(xmax, ae_fabs(x->ptr.p_double[j], _state), _state);
                j = j+jinc;
            }
        }
        *s = *s/tscal;
    }
    
    /*
     * Scale the column norms by 1/TSCAL for return.
     */
    if( ae_fp_neq(tscal,(double)(1)) )
    {
        v = 1/tscal;
        ae_v_muld(&cnorm->ptr.p_double[1], 1, ae_v_len(1,n), v);
    }
}




/*************************************************************************
Real implementation of CMatrixScaledTRSafeSolve

  -- ALGLIB routine --
     21.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixscaledtrsafesolve(/* Real    */ ae_matrix* a,
     double sa,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     ae_bool isupper,
     ae_int_t trans,
     ae_bool isunit,
     double maxgrowth,
     ae_state *_state)
{
    ae_frame _frame_block;
    double lnmax;
    double nrmb;
    double nrmx;
    ae_int_t i;
    ae_complex alpha;
    ae_complex beta;
    double vr;
    ae_complex cx;
    ae_vector tmp;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tmp, 0, DT_REAL, _state);

    ae_assert(n>0, "RMatrixTRSafeSolve: incorrect N!", _state);
    ae_assert(trans==0||trans==1, "RMatrixTRSafeSolve: incorrect Trans!", _state);
    result = ae_true;
    lnmax = ae_log(ae_maxrealnumber, _state);
    
    /*
     * Quick return if possible
     */
    if( n<=0 )
    {
        ae_frame_leave(_state);
        return result;
    }
    
    /*
     * Load norms: right part and X
     */
    nrmb = (double)(0);
    for(i=0; i<=n-1; i++)
    {
        nrmb = ae_maxreal(nrmb, ae_fabs(x->ptr.p_double[i], _state), _state);
    }
    nrmx = (double)(0);
    
    /*
     * Solve
     */
    ae_vector_set_length(&tmp, n, _state);
    result = ae_true;
    if( isupper&&trans==0 )
    {
        
        /*
         * U*x = b
         */
        for(i=n-1; i>=0; i--)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_complex_from_d(a->ptr.pp_double[i][i]*sa);
            }
            if( i<n-1 )
            {
                ae_v_moved(&tmp.ptr.p_double[i+1], 1, &a->ptr.pp_double[i][i+1], 1, ae_v_len(i+1,n-1), sa);
                vr = ae_v_dotproduct(&tmp.ptr.p_double[i+1], 1, &x->ptr.p_double[i+1], 1, ae_v_len(i+1,n-1));
                beta = ae_complex_from_d(x->ptr.p_double[i]-vr);
            }
            else
            {
                beta = ae_complex_from_d(x->ptr.p_double[i]);
            }
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_double[i] = cx.x;
        }
        ae_frame_leave(_state);
        return result;
    }
    if( !isupper&&trans==0 )
    {
        
        /*
         * L*x = b
         */
        for(i=0; i<=n-1; i++)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_complex_from_d(a->ptr.pp_double[i][i]*sa);
            }
            if( i>0 )
            {
                ae_v_moved(&tmp.ptr.p_double[0], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(0,i-1), sa);
                vr = ae_v_dotproduct(&tmp.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,i-1));
                beta = ae_complex_from_d(x->ptr.p_double[i]-vr);
            }
            else
            {
                beta = ae_complex_from_d(x->ptr.p_double[i]);
            }
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_double[i] = cx.x;
        }
        ae_frame_leave(_state);
        return result;
    }
    if( isupper&&trans==1 )
    {
        
        /*
         * U^T*x = b
         */
        for(i=0; i<=n-1; i++)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_complex_from_d(a->ptr.pp_double[i][i]*sa);
            }
            beta = ae_complex_from_d(x->ptr.p_double[i]);
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_double[i] = cx.x;
            
            /*
             * update the rest of right part
             */
            if( i<n-1 )
            {
                vr = cx.x;
                ae_v_moved(&tmp.ptr.p_double[i+1], 1, &a->ptr.pp_double[i][i+1], 1, ae_v_len(i+1,n-1), sa);
                ae_v_subd(&x->ptr.p_double[i+1], 1, &tmp.ptr.p_double[i+1], 1, ae_v_len(i+1,n-1), vr);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    if( !isupper&&trans==1 )
    {
        
        /*
         * L^T*x = b
         */
        for(i=n-1; i>=0; i--)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_complex_from_d(a->ptr.pp_double[i][i]*sa);
            }
            beta = ae_complex_from_d(x->ptr.p_double[i]);
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_double[i] = cx.x;
            
            /*
             * update the rest of right part
             */
            if( i>0 )
            {
                vr = cx.x;
                ae_v_moved(&tmp.ptr.p_double[0], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(0,i-1), sa);
                ae_v_subd(&x->ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,i-1), vr);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    result = ae_false;
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Internal subroutine for safe solution of

    SA*op(A)=b
    
where  A  is  NxN  upper/lower  triangular/unitriangular  matrix, op(A) is
either identity transform, transposition or Hermitian transposition, SA is
a scaling factor such that max(|SA*A[i,j]|) is close to 1.0 in magnutude.

This subroutine  limits  relative  growth  of  solution  (in inf-norm)  by
MaxGrowth,  returning  False  if  growth  exceeds MaxGrowth. Degenerate or
near-degenerate matrices are handled correctly (False is returned) as long
as MaxGrowth is significantly less than MaxRealNumber/norm(b).

  -- ALGLIB routine --
     21.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixscaledtrsafesolve(/* Complex */ ae_matrix* a,
     double sa,
     ae_int_t n,
     /* Complex */ ae_vector* x,
     ae_bool isupper,
     ae_int_t trans,
     ae_bool isunit,
     double maxgrowth,
     ae_state *_state)
{
    ae_frame _frame_block;
    double lnmax;
    double nrmb;
    double nrmx;
    ae_int_t i;
    ae_complex alpha;
    ae_complex beta;
    ae_complex vc;
    ae_vector tmp;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tmp, 0, DT_COMPLEX, _state);

    ae_assert(n>0, "CMatrixTRSafeSolve: incorrect N!", _state);
    ae_assert((trans==0||trans==1)||trans==2, "CMatrixTRSafeSolve: incorrect Trans!", _state);
    result = ae_true;
    lnmax = ae_log(ae_maxrealnumber, _state);
    
    /*
     * Quick return if possible
     */
    if( n<=0 )
    {
        ae_frame_leave(_state);
        return result;
    }
    
    /*
     * Load norms: right part and X
     */
    nrmb = (double)(0);
    for(i=0; i<=n-1; i++)
    {
        nrmb = ae_maxreal(nrmb, ae_c_abs(x->ptr.p_complex[i], _state), _state);
    }
    nrmx = (double)(0);
    
    /*
     * Solve
     */
    ae_vector_set_length(&tmp, n, _state);
    result = ae_true;
    if( isupper&&trans==0 )
    {
        
        /*
         * U*x = b
         */
        for(i=n-1; i>=0; i--)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(a->ptr.pp_complex[i][i],sa);
            }
            if( i<n-1 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[i+1], 1, &a->ptr.pp_complex[i][i+1], 1, "N", ae_v_len(i+1,n-1), sa);
                vc = ae_v_cdotproduct(&tmp.ptr.p_complex[i+1], 1, "N", &x->ptr.p_complex[i+1], 1, "N", ae_v_len(i+1,n-1));
                beta = ae_c_sub(x->ptr.p_complex[i],vc);
            }
            else
            {
                beta = x->ptr.p_complex[i];
            }
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
        }
        ae_frame_leave(_state);
        return result;
    }
    if( !isupper&&trans==0 )
    {
        
        /*
         * L*x = b
         */
        for(i=0; i<=n-1; i++)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(a->ptr.pp_complex[i][i],sa);
            }
            if( i>0 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[0], 1, &a->ptr.pp_complex[i][0], 1, "N", ae_v_len(0,i-1), sa);
                vc = ae_v_cdotproduct(&tmp.ptr.p_complex[0], 1, "N", &x->ptr.p_complex[0], 1, "N", ae_v_len(0,i-1));
                beta = ae_c_sub(x->ptr.p_complex[i],vc);
            }
            else
            {
                beta = x->ptr.p_complex[i];
            }
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
        }
        ae_frame_leave(_state);
        return result;
    }
    if( isupper&&trans==1 )
    {
        
        /*
         * U^T*x = b
         */
        for(i=0; i<=n-1; i++)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(a->ptr.pp_complex[i][i],sa);
            }
            beta = x->ptr.p_complex[i];
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
            
            /*
             * update the rest of right part
             */
            if( i<n-1 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[i+1], 1, &a->ptr.pp_complex[i][i+1], 1, "N", ae_v_len(i+1,n-1), sa);
                ae_v_csubc(&x->ptr.p_complex[i+1], 1, &tmp.ptr.p_complex[i+1], 1, "N", ae_v_len(i+1,n-1), vc);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    if( !isupper&&trans==1 )
    {
        
        /*
         * L^T*x = b
         */
        for(i=n-1; i>=0; i--)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(a->ptr.pp_complex[i][i],sa);
            }
            beta = x->ptr.p_complex[i];
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
            
            /*
             * update the rest of right part
             */
            if( i>0 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[0], 1, &a->ptr.pp_complex[i][0], 1, "N", ae_v_len(0,i-1), sa);
                ae_v_csubc(&x->ptr.p_complex[0], 1, &tmp.ptr.p_complex[0], 1, "N", ae_v_len(0,i-1), vc);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    if( isupper&&trans==2 )
    {
        
        /*
         * U^H*x = b
         */
        for(i=0; i<=n-1; i++)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(ae_c_conj(a->ptr.pp_complex[i][i], _state),sa);
            }
            beta = x->ptr.p_complex[i];
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
            
            /*
             * update the rest of right part
             */
            if( i<n-1 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[i+1], 1, &a->ptr.pp_complex[i][i+1], 1, "Conj", ae_v_len(i+1,n-1), sa);
                ae_v_csubc(&x->ptr.p_complex[i+1], 1, &tmp.ptr.p_complex[i+1], 1, "N", ae_v_len(i+1,n-1), vc);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    if( !isupper&&trans==2 )
    {
        
        /*
         * L^T*x = b
         */
        for(i=n-1; i>=0; i--)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(ae_c_conj(a->ptr.pp_complex[i][i], _state),sa);
            }
            beta = x->ptr.p_complex[i];
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
            
            /*
             * update the rest of right part
             */
            if( i>0 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[0], 1, &a->ptr.pp_complex[i][0], 1, "Conj", ae_v_len(0,i-1), sa);
                ae_v_csubc(&x->ptr.p_complex[0], 1, &tmp.ptr.p_complex[0], 1, "N", ae_v_len(0,i-1), vc);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    result = ae_false;
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
complex basic solver-updater for reduced linear system

    alpha*x[i] = beta

solves this equation and updates it in overlfow-safe manner (keeping track
of relative growth of solution).

Parameters:
    Alpha   -   alpha
    Beta    -   beta
    LnMax   -   precomputed Ln(MaxRealNumber)
    BNorm   -   inf-norm of b (right part of original system)
    MaxGrowth-  maximum growth of norm(x) relative to norm(b)
    XNorm   -   inf-norm of other components of X (which are already processed)
                it is updated by CBasicSolveAndUpdate.
    X       -   solution

  -- ALGLIB routine --
     26.01.2009
     Bochkanov Sergey
*************************************************************************/
static ae_bool safesolve_cbasicsolveandupdate(ae_complex alpha,
     ae_complex beta,
     double lnmax,
     double bnorm,
     double maxgrowth,
     double* xnorm,
     ae_complex* x,
     ae_state *_state)
{
    double v;
    ae_bool result;

    x->x = 0;
    x->y = 0;

    result = ae_false;
    if( ae_c_eq_d(alpha,(double)(0)) )
    {
        return result;
    }
    if( ae_c_neq_d(beta,(double)(0)) )
    {
        
        /*
         * alpha*x[i]=beta
         */
        v = ae_log(ae_c_abs(beta, _state), _state)-ae_log(ae_c_abs(alpha, _state), _state);
        if( ae_fp_greater(v,lnmax) )
        {
            return result;
        }
        *x = ae_c_div(beta,alpha);
    }
    else
    {
        
        /*
         * alpha*x[i]=0
         */
        *x = ae_complex_from_i(0);
    }
    
    /*
     * update NrmX, test growth limit
     */
    *xnorm = ae_maxreal(*xnorm, ae_c_abs(*x, _state), _state);
    if( ae_fp_greater(*xnorm,maxgrowth*bnorm) )
    {
        return result;
    }
    result = ae_true;
    return result;
}




/*************************************************************************
Prepares HPC compuations  of  chunked  gradient with HPCChunkedGradient().
You  have to call this function  before  calling  HPCChunkedGradient() for
a new set of weights. You have to call it only once, see example below:

HOW TO PROCESS DATASET WITH THIS FUNCTION:
    Grad:=0
    HPCPrepareChunkedGradient(Weights, WCount, NTotal, NOut, Buf)
    foreach chunk-of-dataset do
        HPCChunkedGradient(...)
    HPCFinalizeChunkedGradient(Buf, Grad)

*************************************************************************/
void hpcpreparechunkedgradient(/* Real    */ ae_vector* weights,
     ae_int_t wcount,
     ae_int_t ntotal,
     ae_int_t nin,
     ae_int_t nout,
     mlpbuffers* buf,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t batch4size;
    ae_int_t chunksize;


    chunksize = 4;
    batch4size = 3*chunksize*ntotal+chunksize*(2*nout+1);
    if( buf->xy.rows<chunksize||buf->xy.cols<nin+nout )
    {
        ae_matrix_set_length(&buf->xy, chunksize, nin+nout, _state);
    }
    if( buf->xy2.rows<chunksize||buf->xy2.cols<nin+nout )
    {
        ae_matrix_set_length(&buf->xy2, chunksize, nin+nout, _state);
    }
    if( buf->xyrow.cnt<nin+nout )
    {
        ae_vector_set_length(&buf->xyrow, nin+nout, _state);
    }
    if( buf->x.cnt<nin )
    {
        ae_vector_set_length(&buf->x, nin, _state);
    }
    if( buf->y.cnt<nout )
    {
        ae_vector_set_length(&buf->y, nout, _state);
    }
    if( buf->desiredy.cnt<nout )
    {
        ae_vector_set_length(&buf->desiredy, nout, _state);
    }
    if( buf->batch4buf.cnt<batch4size )
    {
        ae_vector_set_length(&buf->batch4buf, batch4size, _state);
    }
    if( buf->hpcbuf.cnt<wcount )
    {
        ae_vector_set_length(&buf->hpcbuf, wcount, _state);
    }
    if( buf->g.cnt<wcount )
    {
        ae_vector_set_length(&buf->g, wcount, _state);
    }
    if( !hpccores_hpcpreparechunkedgradientx(weights, wcount, &buf->hpcbuf, _state) )
    {
        for(i=0; i<=wcount-1; i++)
        {
            buf->hpcbuf.ptr.p_double[i] = 0.0;
        }
    }
    buf->wcount = wcount;
    buf->ntotal = ntotal;
    buf->nin = nin;
    buf->nout = nout;
    buf->chunksize = chunksize;
}


/*************************************************************************
Finalizes HPC compuations  of  chunked gradient with HPCChunkedGradient().
You  have to call this function  after  calling  HPCChunkedGradient()  for
a new set of weights. You have to call it only once, see example below:

HOW TO PROCESS DATASET WITH THIS FUNCTION:
    Grad:=0
    HPCPrepareChunkedGradient(Weights, WCount, NTotal, NOut, Buf)
    foreach chunk-of-dataset do
        HPCChunkedGradient(...)
    HPCFinalizeChunkedGradient(Buf, Grad)

*************************************************************************/
void hpcfinalizechunkedgradient(mlpbuffers* buf,
     /* Real    */ ae_vector* grad,
     ae_state *_state)
{
    ae_int_t i;


    if( !hpccores_hpcfinalizechunkedgradientx(&buf->hpcbuf, buf->wcount, grad, _state) )
    {
        for(i=0; i<=buf->wcount-1; i++)
        {
            grad->ptr.p_double[i] = grad->ptr.p_double[i]+buf->hpcbuf.ptr.p_double[i];
        }
    }
}


/*************************************************************************
Fast kernel for chunked gradient.

*************************************************************************/
ae_bool hpcchunkedgradient(/* Real    */ ae_vector* weights,
     /* Integer */ ae_vector* structinfo,
     /* Real    */ ae_vector* columnmeans,
     /* Real    */ ae_vector* columnsigmas,
     /* Real    */ ae_matrix* xy,
     ae_int_t cstart,
     ae_int_t csize,
     /* Real    */ ae_vector* batch4buf,
     /* Real    */ ae_vector* hpcbuf,
     double* e,
     ae_bool naturalerrorfunc,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_SSE2
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_hpcchunkedgradient(weights, structinfo, columnmeans, columnsigmas, xy, cstart, csize, batch4buf, hpcbuf, e, naturalerrorfunc);
#endif
}


/*************************************************************************
Fast kernel for chunked processing.

*************************************************************************/
ae_bool hpcchunkedprocess(/* Real    */ ae_vector* weights,
     /* Integer */ ae_vector* structinfo,
     /* Real    */ ae_vector* columnmeans,
     /* Real    */ ae_vector* columnsigmas,
     /* Real    */ ae_matrix* xy,
     ae_int_t cstart,
     ae_int_t csize,
     /* Real    */ ae_vector* batch4buf,
     /* Real    */ ae_vector* hpcbuf,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_SSE2
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_hpcchunkedprocess(weights, structinfo, columnmeans, columnsigmas, xy, cstart, csize, batch4buf, hpcbuf);
#endif
}


/*************************************************************************
Stub function.

  -- ALGLIB routine --
     14.06.2013
     Bochkanov Sergey
*************************************************************************/
static ae_bool hpccores_hpcpreparechunkedgradientx(/* Real    */ ae_vector* weights,
     ae_int_t wcount,
     /* Real    */ ae_vector* hpcbuf,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_SSE2
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_hpcpreparechunkedgradientx(weights, wcount, hpcbuf);
#endif
}


/*************************************************************************
Stub function.

  -- ALGLIB routine --
     14.06.2013
     Bochkanov Sergey
*************************************************************************/
static ae_bool hpccores_hpcfinalizechunkedgradientx(/* Real    */ ae_vector* buf,
     ae_int_t wcount,
     /* Real    */ ae_vector* grad,
     ae_state *_state)
{
#ifndef ALGLIB_INTERCEPTS_SSE2
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_hpcfinalizechunkedgradientx(buf, wcount, grad);
#endif
}


void _mlpbuffers_init(void* _p, ae_state *_state)
{
    mlpbuffers *p = (mlpbuffers*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_init(&p->batch4buf, 0, DT_REAL, _state);
    ae_vector_init(&p->hpcbuf, 0, DT_REAL, _state);
    ae_matrix_init(&p->xy, 0, 0, DT_REAL, _state);
    ae_matrix_init(&p->xy2, 0, 0, DT_REAL, _state);
    ae_vector_init(&p->xyrow, 0, DT_REAL, _state);
    ae_vector_init(&p->x, 0, DT_REAL, _state);
    ae_vector_init(&p->y, 0, DT_REAL, _state);
    ae_vector_init(&p->desiredy, 0, DT_REAL, _state);
    ae_vector_init(&p->g, 0, DT_REAL, _state);
    ae_vector_init(&p->tmp0, 0, DT_REAL, _state);
}


void _mlpbuffers_init_copy(void* _dst, void* _src, ae_state *_state)
{
    mlpbuffers *dst = (mlpbuffers*)_dst;
    mlpbuffers *src = (mlpbuffers*)_src;
    dst->chunksize = src->chunksize;
    dst->ntotal = src->ntotal;
    dst->nin = src->nin;
    dst->nout = src->nout;
    dst->wcount = src->wcount;
    ae_vector_init_copy(&dst->batch4buf, &src->batch4buf, _state);
    ae_vector_init_copy(&dst->hpcbuf, &src->hpcbuf, _state);
    ae_matrix_init_copy(&dst->xy, &src->xy, _state);
    ae_matrix_init_copy(&dst->xy2, &src->xy2, _state);
    ae_vector_init_copy(&dst->xyrow, &src->xyrow, _state);
    ae_vector_init_copy(&dst->x, &src->x, _state);
    ae_vector_init_copy(&dst->y, &src->y, _state);
    ae_vector_init_copy(&dst->desiredy, &src->desiredy, _state);
    dst->e = src->e;
    ae_vector_init_copy(&dst->g, &src->g, _state);
    ae_vector_init_copy(&dst->tmp0, &src->tmp0, _state);
}


void _mlpbuffers_clear(void* _p)
{
    mlpbuffers *p = (mlpbuffers*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_clear(&p->batch4buf);
    ae_vector_clear(&p->hpcbuf);
    ae_matrix_clear(&p->xy);
    ae_matrix_clear(&p->xy2);
    ae_vector_clear(&p->xyrow);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->y);
    ae_vector_clear(&p->desiredy);
    ae_vector_clear(&p->g);
    ae_vector_clear(&p->tmp0);
}


void _mlpbuffers_destroy(void* _p)
{
    mlpbuffers *p = (mlpbuffers*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_destroy(&p->batch4buf);
    ae_vector_destroy(&p->hpcbuf);
    ae_matrix_destroy(&p->xy);
    ae_matrix_destroy(&p->xy2);
    ae_vector_destroy(&p->xyrow);
    ae_vector_destroy(&p->x);
    ae_vector_destroy(&p->y);
    ae_vector_destroy(&p->desiredy);
    ae_vector_destroy(&p->g);
    ae_vector_destroy(&p->tmp0);
}




/*************************************************************************
More precise dot-product. Absolute error of  subroutine  result  is  about
1 ulp of max(MX,V), where:
    MX = max( |a[i]*b[i]| )
    V  = |(a,b)|

INPUT PARAMETERS
    A       -   array[0..N-1], vector 1
    B       -   array[0..N-1], vector 2
    N       -   vectors length, N<2^29.
    Temp    -   array[0..N-1], pre-allocated temporary storage

OUTPUT PARAMETERS
    R       -   (A,B)
    RErr    -   estimate of error. This estimate accounts for both  errors
                during  calculation  of  (A,B)  and  errors  introduced by
                rounding of A and B to fit in double (about 1 ulp).

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************/
void xdot(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     ae_int_t n,
     /* Real    */ ae_vector* temp,
     double* r,
     double* rerr,
     ae_state *_state)
{
    ae_int_t i;
    double mx;
    double v;

    *r = 0;
    *rerr = 0;

    
    /*
     * special cases:
     * * N=0
     */
    if( n==0 )
    {
        *r = (double)(0);
        *rerr = (double)(0);
        return;
    }
    mx = (double)(0);
    for(i=0; i<=n-1; i++)
    {
        v = a->ptr.p_double[i]*b->ptr.p_double[i];
        temp->ptr.p_double[i] = v;
        mx = ae_maxreal(mx, ae_fabs(v, _state), _state);
    }
    if( ae_fp_eq(mx,(double)(0)) )
    {
        *r = (double)(0);
        *rerr = (double)(0);
        return;
    }
    xblas_xsum(temp, mx, n, r, rerr, _state);
}


/*************************************************************************
More precise complex dot-product. Absolute error of  subroutine  result is
about 1 ulp of max(MX,V), where:
    MX = max( |a[i]*b[i]| )
    V  = |(a,b)|

INPUT PARAMETERS
    A       -   array[0..N-1], vector 1
    B       -   array[0..N-1], vector 2
    N       -   vectors length, N<2^29.
    Temp    -   array[0..2*N-1], pre-allocated temporary storage

OUTPUT PARAMETERS
    R       -   (A,B)
    RErr    -   estimate of error. This estimate accounts for both  errors
                during  calculation  of  (A,B)  and  errors  introduced by
                rounding of A and B to fit in double (about 1 ulp).

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void xcdot(/* Complex */ ae_vector* a,
     /* Complex */ ae_vector* b,
     ae_int_t n,
     /* Real    */ ae_vector* temp,
     ae_complex* r,
     double* rerr,
     ae_state *_state)
{
    ae_int_t i;
    double mx;
    double v;
    double rerrx;
    double rerry;

    r->x = 0;
    r->y = 0;
    *rerr = 0;

    
    /*
     * special cases:
     * * N=0
     */
    if( n==0 )
    {
        *r = ae_complex_from_i(0);
        *rerr = (double)(0);
        return;
    }
    
    /*
     * calculate real part
     */
    mx = (double)(0);
    for(i=0; i<=n-1; i++)
    {
        v = a->ptr.p_complex[i].x*b->ptr.p_complex[i].x;
        temp->ptr.p_double[2*i+0] = v;
        mx = ae_maxreal(mx, ae_fabs(v, _state), _state);
        v = -a->ptr.p_complex[i].y*b->ptr.p_complex[i].y;
        temp->ptr.p_double[2*i+1] = v;
        mx = ae_maxreal(mx, ae_fabs(v, _state), _state);
    }
    if( ae_fp_eq(mx,(double)(0)) )
    {
        r->x = (double)(0);
        rerrx = (double)(0);
    }
    else
    {
        xblas_xsum(temp, mx, 2*n, &r->x, &rerrx, _state);
    }
    
    /*
     * calculate imaginary part
     */
    mx = (double)(0);
    for(i=0; i<=n-1; i++)
    {
        v = a->ptr.p_complex[i].x*b->ptr.p_complex[i].y;
        temp->ptr.p_double[2*i+0] = v;
        mx = ae_maxreal(mx, ae_fabs(v, _state), _state);
        v = a->ptr.p_complex[i].y*b->ptr.p_complex[i].x;
        temp->ptr.p_double[2*i+1] = v;
        mx = ae_maxreal(mx, ae_fabs(v, _state), _state);
    }
    if( ae_fp_eq(mx,(double)(0)) )
    {
        r->y = (double)(0);
        rerry = (double)(0);
    }
    else
    {
        xblas_xsum(temp, mx, 2*n, &r->y, &rerry, _state);
    }
    
    /*
     * total error
     */
    if( ae_fp_eq(rerrx,(double)(0))&&ae_fp_eq(rerry,(double)(0)) )
    {
        *rerr = (double)(0);
    }
    else
    {
        *rerr = ae_maxreal(rerrx, rerry, _state)*ae_sqrt(1+ae_sqr(ae_minreal(rerrx, rerry, _state)/ae_maxreal(rerrx, rerry, _state), _state), _state);
    }
}


/*************************************************************************
Internal subroutine for extra-precise calculation of SUM(w[i]).

INPUT PARAMETERS:
    W   -   array[0..N-1], values to be added
            W is modified during calculations.
    MX  -   max(W[i])
    N   -   array size
    
OUTPUT PARAMETERS:
    R   -   SUM(w[i])
    RErr-   error estimate for R

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************/
static void xblas_xsum(/* Real    */ ae_vector* w,
     double mx,
     ae_int_t n,
     double* r,
     double* rerr,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;
    ae_int_t ks;
    double v;
    double s;
    double ln2;
    double chunk;
    double invchunk;
    ae_bool allzeros;

    *r = 0;
    *rerr = 0;

    
    /*
     * special cases:
     * * N=0
     * * N is too large to use integer arithmetics
     */
    if( n==0 )
    {
        *r = (double)(0);
        *rerr = (double)(0);
        return;
    }
    if( ae_fp_eq(mx,(double)(0)) )
    {
        *r = (double)(0);
        *rerr = (double)(0);
        return;
    }
    ae_assert(n<536870912, "XDot: N is too large!", _state);
    
    /*
     * Prepare
     */
    ln2 = ae_log((double)(2), _state);
    *rerr = mx*ae_machineepsilon;
    
    /*
     * 1. find S such that 0.5<=S*MX<1
     * 2. multiply W by S, so task is normalized in some sense
     * 3. S:=1/S so we can obtain original vector multiplying by S
     */
    k = ae_round(ae_log(mx, _state)/ln2, _state);
    s = xblas_xfastpow((double)(2), -k, _state);
    while(ae_fp_greater_eq(s*mx,(double)(1)))
    {
        s = 0.5*s;
    }
    while(ae_fp_less(s*mx,0.5))
    {
        s = 2*s;
    }
    ae_v_muld(&w->ptr.p_double[0], 1, ae_v_len(0,n-1), s);
    s = 1/s;
    
    /*
     * find Chunk=2^M such that N*Chunk<2^29
     *
     * we have chosen upper limit (2^29) with enough space left
     * to tolerate possible problems with rounding and N's close
     * to the limit, so we don't want to be very strict here.
     */
    k = ae_trunc(ae_log((double)536870912/(double)n, _state)/ln2, _state);
    chunk = xblas_xfastpow((double)(2), k, _state);
    if( ae_fp_less(chunk,(double)(2)) )
    {
        chunk = (double)(2);
    }
    invchunk = 1/chunk;
    
    /*
     * calculate result
     */
    *r = (double)(0);
    ae_v_muld(&w->ptr.p_double[0], 1, ae_v_len(0,n-1), chunk);
    for(;;)
    {
        s = s*invchunk;
        allzeros = ae_true;
        ks = 0;
        for(i=0; i<=n-1; i++)
        {
            v = w->ptr.p_double[i];
            k = ae_trunc(v, _state);
            if( ae_fp_neq(v,(double)(k)) )
            {
                allzeros = ae_false;
            }
            w->ptr.p_double[i] = chunk*(v-k);
            ks = ks+k;
        }
        *r = *r+s*ks;
        v = ae_fabs(*r, _state);
        if( allzeros||ae_fp_eq(s*n+mx,mx) )
        {
            break;
        }
    }
    
    /*
     * correct error
     */
    *rerr = ae_maxreal(*rerr, ae_fabs(*r, _state)*ae_machineepsilon, _state);
}


/*************************************************************************
Fast Pow

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************/
static double xblas_xfastpow(double r, ae_int_t n, ae_state *_state)
{
    double result;


    result = (double)(0);
    if( n>0 )
    {
        if( n%2==0 )
        {
            result = ae_sqr(xblas_xfastpow(r, n/2, _state), _state);
        }
        else
        {
            result = r*xblas_xfastpow(r, n-1, _state);
        }
        return result;
    }
    if( n==0 )
    {
        result = (double)(1);
    }
    if( n<0 )
    {
        result = xblas_xfastpow(1/r, -n, _state);
    }
    return result;
}




/*************************************************************************
Normalizes direction/step pair: makes |D|=1, scales Stp.
If |D|=0, it returns, leavind D/Stp unchanged.

  -- ALGLIB --
     Copyright 01.04.2010 by Bochkanov Sergey
*************************************************************************/
void linminnormalized(/* Real    */ ae_vector* d,
     double* stp,
     ae_int_t n,
     ae_state *_state)
{
    double mx;
    double s;
    ae_int_t i;


    
    /*
     * first, scale D to avoid underflow/overflow durng squaring
     */
    mx = (double)(0);
    for(i=0; i<=n-1; i++)
    {
        mx = ae_maxreal(mx, ae_fabs(d->ptr.p_double[i], _state), _state);
    }
    if( ae_fp_eq(mx,(double)(0)) )
    {
        return;
    }
    s = 1/mx;
    ae_v_muld(&d->ptr.p_double[0], 1, ae_v_len(0,n-1), s);
    *stp = *stp/s;
    
    /*
     * normalize D
     */
    s = ae_v_dotproduct(&d->ptr.p_double[0], 1, &d->ptr.p_double[0], 1, ae_v_len(0,n-1));
    s = 1/ae_sqrt(s, _state);
    ae_v_muld(&d->ptr.p_double[0], 1, ae_v_len(0,n-1), s);
    *stp = *stp/s;
}


/*************************************************************************
THE  PURPOSE  OF  MCSRCH  IS  TO  FIND A STEP WHICH SATISFIES A SUFFICIENT
DECREASE CONDITION AND A CURVATURE CONDITION.

AT EACH STAGE THE SUBROUTINE  UPDATES  AN  INTERVAL  OF  UNCERTAINTY  WITH
ENDPOINTS  STX  AND  STY.  THE INTERVAL OF UNCERTAINTY IS INITIALLY CHOSEN
SO THAT IT CONTAINS A MINIMIZER OF THE MODIFIED FUNCTION

    F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).

IF  A STEP  IS OBTAINED FOR  WHICH THE MODIFIED FUNCTION HAS A NONPOSITIVE
FUNCTION  VALUE  AND  NONNEGATIVE  DERIVATIVE,   THEN   THE   INTERVAL  OF
UNCERTAINTY IS CHOSEN SO THAT IT CONTAINS A MINIMIZER OF F(X+STP*S).

THE  ALGORITHM  IS  DESIGNED TO FIND A STEP WHICH SATISFIES THE SUFFICIENT
DECREASE CONDITION

    F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),

AND THE CURVATURE CONDITION

    ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).

IF  FTOL  IS  LESS  THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION IS BOUNDED
BELOW,  THEN  THERE  IS  ALWAYS  A  STEP  WHICH SATISFIES BOTH CONDITIONS.
IF  NO  STEP  CAN BE FOUND  WHICH  SATISFIES  BOTH  CONDITIONS,  THEN  THE
ALGORITHM  USUALLY STOPS  WHEN  ROUNDING ERRORS  PREVENT FURTHER PROGRESS.
IN THIS CASE STP ONLY SATISFIES THE SUFFICIENT DECREASE CONDITION.


:::::::::::::IMPORTANT NOTES:::::::::::::

NOTE 1:

This routine  guarantees that it will stop at the last point where function
value was calculated. It won't make several additional function evaluations
after finding good point. So if you store function evaluations requested by
this routine, you can be sure that last one is the point where we've stopped.

NOTE 2:

when 0<StpMax<StpMin, algorithm will terminate with INFO=5 and Stp=StpMax

NOTE 3:

this algorithm guarantees that, if MCINFO=1 or MCINFO=5, then:
* F(final_point)<F(initial_point) - strict inequality
* final_point<>initial_point - after rounding to machine precision
:::::::::::::::::::::::::::::::::::::::::


PARAMETERS DESCRIPRION

STAGE IS ZERO ON FIRST CALL, ZERO ON FINAL EXIT

N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF VARIABLES.

X IS  AN  ARRAY  OF  LENGTH N. ON INPUT IT MUST CONTAIN THE BASE POINT FOR
THE LINE SEARCH. ON OUTPUT IT CONTAINS X+STP*S.

F IS  A  VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F AT X. ON OUTPUT
IT CONTAINS THE VALUE OF F AT X + STP*S.

G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE GRADIENT OF F AT X.
ON OUTPUT IT CONTAINS THE GRADIENT OF F AT X + STP*S.

S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE SEARCH DIRECTION.

STP  IS  A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN INITIAL ESTIMATE
OF A SATISFACTORY STEP. ON OUTPUT STP CONTAINS THE FINAL ESTIMATE.

FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION OCCURS WHEN THE
SUFFICIENT DECREASE CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
SATISFIED.

XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS WHEN THE RELATIVE
WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL.

STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH SPECIFY LOWER  AND
UPPER BOUNDS FOR THE STEP.

MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION OCCURS WHEN THE
NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN ITERATION.

INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
    INFO = 0  IMPROPER INPUT PARAMETERS.

    INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
              DIRECTIONAL DERIVATIVE CONDITION HOLD.

    INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
              IS AT MOST XTOL.

    INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.

    INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.

    INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.

    INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
              THERE MAY NOT BE A STEP WHICH SATISFIES THE
              SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
              TOLERANCES MAY BE TOO SMALL.

NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.

WA IS A WORK ARRAY OF LENGTH N.

ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
JORGE J. MORE', DAVID J. THUENTE
*************************************************************************/
void mcsrch(ae_int_t n,
     /* Real    */ ae_vector* x,
     double* f,
     /* Real    */ ae_vector* g,
     /* Real    */ ae_vector* s,
     double* stp,
     double stpmax,
     double gtol,
     ae_int_t* info,
     ae_int_t* nfev,
     /* Real    */ ae_vector* wa,
     linminstate* state,
     ae_int_t* stage,
     ae_state *_state)
{
    ae_int_t i;
    double v;
    double p5;
    double p66;
    double zero;


    
    /*
     * init
     */
    p5 = 0.5;
    p66 = 0.66;
    state->xtrapf = 4.0;
    zero = (double)(0);
    if( ae_fp_eq(stpmax,(double)(0)) )
    {
        stpmax = linmin_defstpmax;
    }
    if( ae_fp_less(*stp,linmin_stpmin) )
    {
        *stp = linmin_stpmin;
    }
    if( ae_fp_greater(*stp,stpmax) )
    {
        *stp = stpmax;
    }
    
    /*
     * Main cycle
     */
    for(;;)
    {
        if( *stage==0 )
        {
            
            /*
             * NEXT
             */
            *stage = 2;
            continue;
        }
        if( *stage==2 )
        {
            state->infoc = 1;
            *info = 0;
            
            /*
             *     CHECK THE INPUT PARAMETERS FOR ERRORS.
             */
            if( ae_fp_less(stpmax,linmin_stpmin)&&ae_fp_greater(stpmax,(double)(0)) )
            {
                *info = 5;
                *stp = stpmax;
                *stage = 0;
                return;
            }
            if( ((((((n<=0||ae_fp_less_eq(*stp,(double)(0)))||ae_fp_less(linmin_ftol,(double)(0)))||ae_fp_less(gtol,zero))||ae_fp_less(linmin_xtol,zero))||ae_fp_less(linmin_stpmin,zero))||ae_fp_less(stpmax,linmin_stpmin))||linmin_maxfev<=0 )
            {
                *stage = 0;
                return;
            }
            
            /*
             *     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
             *     AND CHECK THAT S IS A DESCENT DIRECTION.
             */
            v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &s->ptr.p_double[0], 1, ae_v_len(0,n-1));
            state->dginit = v;
            if( ae_fp_greater_eq(state->dginit,(double)(0)) )
            {
                *stage = 0;
                return;
            }
            
            /*
             *     INITIALIZE LOCAL VARIABLES.
             */
            state->brackt = ae_false;
            state->stage1 = ae_true;
            *nfev = 0;
            state->finit = *f;
            state->dgtest = linmin_ftol*state->dginit;
            state->width = stpmax-linmin_stpmin;
            state->width1 = state->width/p5;
            ae_v_move(&wa->ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
            
            /*
             *     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
             *     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
             *     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
             *     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
             *     THE INTERVAL OF UNCERTAINTY.
             *     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
             *     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
             */
            state->stx = (double)(0);
            state->fx = state->finit;
            state->dgx = state->dginit;
            state->sty = (double)(0);
            state->fy = state->finit;
            state->dgy = state->dginit;
            
            /*
             * NEXT
             */
            *stage = 3;
            continue;
        }
        if( *stage==3 )
        {
            
            /*
             *     START OF ITERATION.
             *
             *     SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
             *     TO THE PRESENT INTERVAL OF UNCERTAINTY.
             */
            if( state->brackt )
            {
                if( ae_fp_less(state->stx,state->sty) )
                {
                    state->stmin = state->stx;
                    state->stmax = state->sty;
                }
                else
                {
                    state->stmin = state->sty;
                    state->stmax = state->stx;
                }
            }
            else
            {
                state->stmin = state->stx;
                state->stmax = *stp+state->xtrapf*(*stp-state->stx);
            }
            
            /*
             *        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
             */
            if( ae_fp_greater(*stp,stpmax) )
            {
                *stp = stpmax;
            }
            if( ae_fp_less(*stp,linmin_stpmin) )
            {
                *stp = linmin_stpmin;
            }
            
            /*
             *        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
             *        STP BE THE LOWEST POINT OBTAINED SO FAR.
             */
            if( (((state->brackt&&(ae_fp_less_eq(*stp,state->stmin)||ae_fp_greater_eq(*stp,state->stmax)))||*nfev>=linmin_maxfev-1)||state->infoc==0)||(state->brackt&&ae_fp_less_eq(state->stmax-state->stmin,linmin_xtol*state->stmax)) )
            {
                *stp = state->stx;
            }
            
            /*
             *        EVALUATE THE FUNCTION AND GRADIENT AT STP
             *        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
             */
            ae_v_move(&x->ptr.p_double[0], 1, &wa->ptr.p_double[0], 1, ae_v_len(0,n-1));
            ae_v_addd(&x->ptr.p_double[0], 1, &s->ptr.p_double[0], 1, ae_v_len(0,n-1), *stp);
            
            /*
             * NEXT
             */
            *stage = 4;
            return;
        }
        if( *stage==4 )
        {
            *info = 0;
            *nfev = *nfev+1;
            v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &s->ptr.p_double[0], 1, ae_v_len(0,n-1));
            state->dg = v;
            state->ftest1 = state->finit+*stp*state->dgtest;
            
            /*
             *        TEST FOR CONVERGENCE.
             */
            if( (state->brackt&&(ae_fp_less_eq(*stp,state->stmin)||ae_fp_greater_eq(*stp,state->stmax)))||state->infoc==0 )
            {
                *info = 6;
            }
            if( ((ae_fp_eq(*stp,stpmax)&&ae_fp_less(*f,state->finit))&&ae_fp_less_eq(*f,state->ftest1))&&ae_fp_less_eq(state->dg,state->dgtest) )
            {
                *info = 5;
            }
            if( ae_fp_eq(*stp,linmin_stpmin)&&((ae_fp_greater_eq(*f,state->finit)||ae_fp_greater(*f,state->ftest1))||ae_fp_greater_eq(state->dg,state->dgtest)) )
            {
                *info = 4;
            }
            if( *nfev>=linmin_maxfev )
            {
                *info = 3;
            }
            if( state->brackt&&ae_fp_less_eq(state->stmax-state->stmin,linmin_xtol*state->stmax) )
            {
                *info = 2;
            }
            if( (ae_fp_less(*f,state->finit)&&ae_fp_less_eq(*f,state->ftest1))&&ae_fp_less_eq(ae_fabs(state->dg, _state),-gtol*state->dginit) )
            {
                *info = 1;
            }
            
            /*
             *        CHECK FOR TERMINATION.
             */
            if( *info!=0 )
            {
                
                /*
                 * Check guarantees provided by the function for INFO=1 or INFO=5
                 */
                if( *info==1||*info==5 )
                {
                    v = 0.0;
                    for(i=0; i<=n-1; i++)
                    {
                        v = v+(wa->ptr.p_double[i]-x->ptr.p_double[i])*(wa->ptr.p_double[i]-x->ptr.p_double[i]);
                    }
                    if( ae_fp_greater_eq(*f,state->finit)||ae_fp_eq(v,0.0) )
                    {
                        *info = 6;
                    }
                }
                *stage = 0;
                return;
            }
            
            /*
             *        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
             *        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
             */
            if( (state->stage1&&ae_fp_less_eq(*f,state->ftest1))&&ae_fp_greater_eq(state->dg,ae_minreal(linmin_ftol, gtol, _state)*state->dginit) )
            {
                state->stage1 = ae_false;
            }
            
            /*
             *        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
             *        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
             *        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
             *        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
             *        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
             */
            if( (state->stage1&&ae_fp_less_eq(*f,state->fx))&&ae_fp_greater(*f,state->ftest1) )
            {
                
                /*
                 *           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
                 */
                state->fm = *f-*stp*state->dgtest;
                state->fxm = state->fx-state->stx*state->dgtest;
                state->fym = state->fy-state->sty*state->dgtest;
                state->dgm = state->dg-state->dgtest;
                state->dgxm = state->dgx-state->dgtest;
                state->dgym = state->dgy-state->dgtest;
                
                /*
                 *           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                 *           AND TO COMPUTE THE NEW STEP.
                 */
                linmin_mcstep(&state->stx, &state->fxm, &state->dgxm, &state->sty, &state->fym, &state->dgym, stp, state->fm, state->dgm, &state->brackt, state->stmin, state->stmax, &state->infoc, _state);
                
                /*
                 *           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
                 */
                state->fx = state->fxm+state->stx*state->dgtest;
                state->fy = state->fym+state->sty*state->dgtest;
                state->dgx = state->dgxm+state->dgtest;
                state->dgy = state->dgym+state->dgtest;
            }
            else
            {
                
                /*
                 *           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                 *           AND TO COMPUTE THE NEW STEP.
                 */
                linmin_mcstep(&state->stx, &state->fx, &state->dgx, &state->sty, &state->fy, &state->dgy, stp, *f, state->dg, &state->brackt, state->stmin, state->stmax, &state->infoc, _state);
            }
            
            /*
             *        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
             *        INTERVAL OF UNCERTAINTY.
             */
            if( state->brackt )
            {
                if( ae_fp_greater_eq(ae_fabs(state->sty-state->stx, _state),p66*state->width1) )
                {
                    *stp = state->stx+p5*(state->sty-state->stx);
                }
                state->width1 = state->width;
                state->width = ae_fabs(state->sty-state->stx, _state);
            }
            
            /*
             *  NEXT.
             */
            *stage = 3;
            continue;
        }
    }
}


/*************************************************************************
These functions perform Armijo line search using  at  most  FMAX  function
evaluations.  It  doesn't  enforce  some  kind  of  " sufficient decrease"
criterion - it just tries different Armijo steps and returns optimum found
so far.

Optimization is done using F-rcomm interface:
* ArmijoCreate initializes State structure
  (reusing previously allocated buffers)
* ArmijoIteration is subsequently called
* ArmijoResults returns results

INPUT PARAMETERS:
    N       -   problem size
    X       -   array[N], starting point
    F       -   F(X+S*STP)
    S       -   step direction, S>0
    STP     -   step length
    STPMAX  -   maximum value for STP or zero (if no limit is imposed)
    FMAX    -   maximum number of function evaluations
    State   -   optimization state

  -- ALGLIB --
     Copyright 05.10.2010 by Bochkanov Sergey
*************************************************************************/
void armijocreate(ae_int_t n,
     /* Real    */ ae_vector* x,
     double f,
     /* Real    */ ae_vector* s,
     double stp,
     double stpmax,
     ae_int_t fmax,
     armijostate* state,
     ae_state *_state)
{


    if( state->x.cnt<n )
    {
        ae_vector_set_length(&state->x, n, _state);
    }
    if( state->xbase.cnt<n )
    {
        ae_vector_set_length(&state->xbase, n, _state);
    }
    if( state->s.cnt<n )
    {
        ae_vector_set_length(&state->s, n, _state);
    }
    state->stpmax = stpmax;
    state->fmax = fmax;
    state->stplen = stp;
    state->fcur = f;
    state->n = n;
    ae_v_move(&state->xbase.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->s.ptr.p_double[0], 1, &s->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_vector_set_length(&state->rstate.ia, 0+1, _state);
    ae_vector_set_length(&state->rstate.ra, 0+1, _state);
    state->rstate.stage = -1;
}


/*************************************************************************
This is rcomm-based search function

  -- ALGLIB --
     Copyright 05.10.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool armijoiteration(armijostate* state, ae_state *_state)
{
    double v;
    ae_int_t n;
    ae_bool result;


    
    /*
     * Reverse communication preparations
     * I know it looks ugly, but it works the same way
     * anywhere from C++ to Python.
     *
     * This code initializes locals by:
     * * random values determined during code
     *   generation - on first subroutine call
     * * values from previous call - on subsequent calls
     */
    if( state->rstate.stage>=0 )
    {
        n = state->rstate.ia.ptr.p_int[0];
        v = state->rstate.ra.ptr.p_double[0];
    }
    else
    {
        n = -983;
        v = -989;
    }
    if( state->rstate.stage==0 )
    {
        goto lbl_0;
    }
    if( state->rstate.stage==1 )
    {
        goto lbl_1;
    }
    if( state->rstate.stage==2 )
    {
        goto lbl_2;
    }
    if( state->rstate.stage==3 )
    {
        goto lbl_3;
    }
    
    /*
     * Routine body
     */
    if( (ae_fp_less_eq(state->stplen,(double)(0))||ae_fp_less(state->stpmax,(double)(0)))||state->fmax<2 )
    {
        state->info = 0;
        result = ae_false;
        return result;
    }
    if( ae_fp_less_eq(state->stplen,linmin_stpmin) )
    {
        state->info = 4;
        result = ae_false;
        return result;
    }
    n = state->n;
    state->nfev = 0;
    
    /*
     * We always need F
     */
    state->needf = ae_true;
    
    /*
     * Bound StpLen
     */
    if( ae_fp_greater(state->stplen,state->stpmax)&&ae_fp_neq(state->stpmax,(double)(0)) )
    {
        state->stplen = state->stpmax;
    }
    
    /*
     * Increase length
     */
    v = state->stplen*linmin_armijofactor;
    if( ae_fp_greater(v,state->stpmax)&&ae_fp_neq(state->stpmax,(double)(0)) )
    {
        v = state->stpmax;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->x.ptr.p_double[0], 1, &state->s.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->nfev = state->nfev+1;
    if( ae_fp_greater_eq(state->f,state->fcur) )
    {
        goto lbl_4;
    }
    state->stplen = v;
    state->fcur = state->f;
lbl_6:
    if( ae_false )
    {
        goto lbl_7;
    }
    
    /*
     * test stopping conditions
     */
    if( state->nfev>=state->fmax )
    {
        state->info = 3;
        result = ae_false;
        return result;
    }
    if( ae_fp_greater_eq(state->stplen,state->stpmax) )
    {
        state->info = 5;
        result = ae_false;
        return result;
    }
    
    /*
     * evaluate F
     */
    v = state->stplen*linmin_armijofactor;
    if( ae_fp_greater(v,state->stpmax)&&ae_fp_neq(state->stpmax,(double)(0)) )
    {
        v = state->stpmax;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->x.ptr.p_double[0], 1, &state->s.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    state->nfev = state->nfev+1;
    
    /*
     * make decision
     */
    if( ae_fp_less(state->f,state->fcur) )
    {
        state->stplen = v;
        state->fcur = state->f;
    }
    else
    {
        state->info = 1;
        result = ae_false;
        return result;
    }
    goto lbl_6;
lbl_7:
lbl_4:
    
    /*
     * Decrease length
     */
    v = state->stplen/linmin_armijofactor;
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->x.ptr.p_double[0], 1, &state->s.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    state->nfev = state->nfev+1;
    if( ae_fp_greater_eq(state->f,state->fcur) )
    {
        goto lbl_8;
    }
    state->stplen = state->stplen/linmin_armijofactor;
    state->fcur = state->f;
lbl_10:
    if( ae_false )
    {
        goto lbl_11;
    }
    
    /*
     * test stopping conditions
     */
    if( state->nfev>=state->fmax )
    {
        state->info = 3;
        result = ae_false;
        return result;
    }
    if( ae_fp_less_eq(state->stplen,linmin_stpmin) )
    {
        state->info = 4;
        result = ae_false;
        return result;
    }
    
    /*
     * evaluate F
     */
    v = state->stplen/linmin_armijofactor;
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->x.ptr.p_double[0], 1, &state->s.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->nfev = state->nfev+1;
    
    /*
     * make decision
     */
    if( ae_fp_less(state->f,state->fcur) )
    {
        state->stplen = state->stplen/linmin_armijofactor;
        state->fcur = state->f;
    }
    else
    {
        state->info = 1;
        result = ae_false;
        return result;
    }
    goto lbl_10;
lbl_11:
lbl_8:
    
    /*
     * Nothing to be done
     */
    state->info = 1;
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ra.ptr.p_double[0] = v;
    return result;
}


/*************************************************************************
Results of Armijo search

OUTPUT PARAMETERS:
    INFO    -   on output it is set to one of the return codes:
                * 0     improper input params
                * 1     optimum step is found with at most FMAX evaluations
                * 3     FMAX evaluations were used,
                        X contains optimum found so far
                * 4     step is at lower bound STPMIN
                * 5     step is at upper bound
    STP     -   step length (in case of failure it is still returned)
    F       -   function value (in case of failure it is still returned)

  -- ALGLIB --
     Copyright 05.10.2010 by Bochkanov Sergey
*************************************************************************/
void armijoresults(armijostate* state,
     ae_int_t* info,
     double* stp,
     double* f,
     ae_state *_state)
{


    *info = state->info;
    *stp = state->stplen;
    *f = state->fcur;
}


static void linmin_mcstep(double* stx,
     double* fx,
     double* dx,
     double* sty,
     double* fy,
     double* dy,
     double* stp,
     double fp,
     double dp,
     ae_bool* brackt,
     double stmin,
     double stmax,
     ae_int_t* info,
     ae_state *_state)
{
    ae_bool bound;
    double gamma;
    double p;
    double q;
    double r;
    double s;
    double sgnd;
    double stpc;
    double stpf;
    double stpq;
    double theta;


    *info = 0;
    
    /*
     *     CHECK THE INPUT PARAMETERS FOR ERRORS.
     */
    if( ((*brackt&&(ae_fp_less_eq(*stp,ae_minreal(*stx, *sty, _state))||ae_fp_greater_eq(*stp,ae_maxreal(*stx, *sty, _state))))||ae_fp_greater_eq(*dx*(*stp-(*stx)),(double)(0)))||ae_fp_less(stmax,stmin) )
    {
        return;
    }
    
    /*
     *     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
     */
    sgnd = dp*(*dx/ae_fabs(*dx, _state));
    
    /*
     *     FIRST CASE. A HIGHER FUNCTION VALUE.
     *     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
     *     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
     *     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
     */
    if( ae_fp_greater(fp,*fx) )
    {
        *info = 1;
        bound = ae_true;
        theta = 3*(*fx-fp)/(*stp-(*stx))+(*dx)+dp;
        s = ae_maxreal(ae_fabs(theta, _state), ae_maxreal(ae_fabs(*dx, _state), ae_fabs(dp, _state), _state), _state);
        gamma = s*ae_sqrt(ae_sqr(theta/s, _state)-*dx/s*(dp/s), _state);
        if( ae_fp_less(*stp,*stx) )
        {
            gamma = -gamma;
        }
        p = gamma-(*dx)+theta;
        q = gamma-(*dx)+gamma+dp;
        r = p/q;
        stpc = *stx+r*(*stp-(*stx));
        stpq = *stx+*dx/((*fx-fp)/(*stp-(*stx))+(*dx))/2*(*stp-(*stx));
        if( ae_fp_less(ae_fabs(stpc-(*stx), _state),ae_fabs(stpq-(*stx), _state)) )
        {
            stpf = stpc;
        }
        else
        {
            stpf = stpc+(stpq-stpc)/2;
        }
        *brackt = ae_true;
    }
    else
    {
        if( ae_fp_less(sgnd,(double)(0)) )
        {
            
            /*
             *     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
             *     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
             *     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
             *     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
             */
            *info = 2;
            bound = ae_false;
            theta = 3*(*fx-fp)/(*stp-(*stx))+(*dx)+dp;
            s = ae_maxreal(ae_fabs(theta, _state), ae_maxreal(ae_fabs(*dx, _state), ae_fabs(dp, _state), _state), _state);
            gamma = s*ae_sqrt(ae_sqr(theta/s, _state)-*dx/s*(dp/s), _state);
            if( ae_fp_greater(*stp,*stx) )
            {
                gamma = -gamma;
            }
            p = gamma-dp+theta;
            q = gamma-dp+gamma+(*dx);
            r = p/q;
            stpc = *stp+r*(*stx-(*stp));
            stpq = *stp+dp/(dp-(*dx))*(*stx-(*stp));
            if( ae_fp_greater(ae_fabs(stpc-(*stp), _state),ae_fabs(stpq-(*stp), _state)) )
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
            *brackt = ae_true;
        }
        else
        {
            if( ae_fp_less(ae_fabs(dp, _state),ae_fabs(*dx, _state)) )
            {
                
                /*
                 *     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                 *     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
                 *     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
                 *     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
                 *     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
                 *     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
                 *     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
                 *     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
                 */
                *info = 3;
                bound = ae_true;
                theta = 3*(*fx-fp)/(*stp-(*stx))+(*dx)+dp;
                s = ae_maxreal(ae_fabs(theta, _state), ae_maxreal(ae_fabs(*dx, _state), ae_fabs(dp, _state), _state), _state);
                
                /*
                 *        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
                 *        TO INFINITY IN THE DIRECTION OF THE STEP.
                 */
                gamma = s*ae_sqrt(ae_maxreal((double)(0), ae_sqr(theta/s, _state)-*dx/s*(dp/s), _state), _state);
                if( ae_fp_greater(*stp,*stx) )
                {
                    gamma = -gamma;
                }
                p = gamma-dp+theta;
                q = gamma+(*dx-dp)+gamma;
                r = p/q;
                if( ae_fp_less(r,(double)(0))&&ae_fp_neq(gamma,(double)(0)) )
                {
                    stpc = *stp+r*(*stx-(*stp));
                }
                else
                {
                    if( ae_fp_greater(*stp,*stx) )
                    {
                        stpc = stmax;
                    }
                    else
                    {
                        stpc = stmin;
                    }
                }
                stpq = *stp+dp/(dp-(*dx))*(*stx-(*stp));
                if( *brackt )
                {
                    if( ae_fp_less(ae_fabs(*stp-stpc, _state),ae_fabs(*stp-stpq, _state)) )
                    {
                        stpf = stpc;
                    }
                    else
                    {
                        stpf = stpq;
                    }
                }
                else
                {
                    if( ae_fp_greater(ae_fabs(*stp-stpc, _state),ae_fabs(*stp-stpq, _state)) )
                    {
                        stpf = stpc;
                    }
                    else
                    {
                        stpf = stpq;
                    }
                }
            }
            else
            {
                
                /*
                 *     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                 *     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
                 *     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
                 *     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
                 */
                *info = 4;
                bound = ae_false;
                if( *brackt )
                {
                    theta = 3*(fp-(*fy))/(*sty-(*stp))+(*dy)+dp;
                    s = ae_maxreal(ae_fabs(theta, _state), ae_maxreal(ae_fabs(*dy, _state), ae_fabs(dp, _state), _state), _state);
                    gamma = s*ae_sqrt(ae_sqr(theta/s, _state)-*dy/s*(dp/s), _state);
                    if( ae_fp_greater(*stp,*sty) )
                    {
                        gamma = -gamma;
                    }
                    p = gamma-dp+theta;
                    q = gamma-dp+gamma+(*dy);
                    r = p/q;
                    stpc = *stp+r*(*sty-(*stp));
                    stpf = stpc;
                }
                else
                {
                    if( ae_fp_greater(*stp,*stx) )
                    {
                        stpf = stmax;
                    }
                    else
                    {
                        stpf = stmin;
                    }
                }
            }
        }
    }
    
    /*
     *     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
     *     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
     */
    if( ae_fp_greater(fp,*fx) )
    {
        *sty = *stp;
        *fy = fp;
        *dy = dp;
    }
    else
    {
        if( ae_fp_less(sgnd,0.0) )
        {
            *sty = *stx;
            *fy = *fx;
            *dy = *dx;
        }
        *stx = *stp;
        *fx = fp;
        *dx = dp;
    }
    
    /*
     *     COMPUTE THE NEW STEP AND SAFEGUARD IT.
     */
    stpf = ae_minreal(stmax, stpf, _state);
    stpf = ae_maxreal(stmin, stpf, _state);
    *stp = stpf;
    if( *brackt&&bound )
    {
        if( ae_fp_greater(*sty,*stx) )
        {
            *stp = ae_minreal(*stx+0.66*(*sty-(*stx)), *stp, _state);
        }
        else
        {
            *stp = ae_maxreal(*stx+0.66*(*sty-(*stx)), *stp, _state);
        }
    }
}


void _linminstate_init(void* _p, ae_state *_state)
{
    linminstate *p = (linminstate*)_p;
    ae_touch_ptr((void*)p);
}


void _linminstate_init_copy(void* _dst, void* _src, ae_state *_state)
{
    linminstate *dst = (linminstate*)_dst;
    linminstate *src = (linminstate*)_src;
    dst->brackt = src->brackt;
    dst->stage1 = src->stage1;
    dst->infoc = src->infoc;
    dst->dg = src->dg;
    dst->dgm = src->dgm;
    dst->dginit = src->dginit;
    dst->dgtest = src->dgtest;
    dst->dgx = src->dgx;
    dst->dgxm = src->dgxm;
    dst->dgy = src->dgy;
    dst->dgym = src->dgym;
    dst->finit = src->finit;
    dst->ftest1 = src->ftest1;
    dst->fm = src->fm;
    dst->fx = src->fx;
    dst->fxm = src->fxm;
    dst->fy = src->fy;
    dst->fym = src->fym;
    dst->stx = src->stx;
    dst->sty = src->sty;
    dst->stmin = src->stmin;
    dst->stmax = src->stmax;
    dst->width = src->width;
    dst->width1 = src->width1;
    dst->xtrapf = src->xtrapf;
}


void _linminstate_clear(void* _p)
{
    linminstate *p = (linminstate*)_p;
    ae_touch_ptr((void*)p);
}


void _linminstate_destroy(void* _p)
{
    linminstate *p = (linminstate*)_p;
    ae_touch_ptr((void*)p);
}


void _armijostate_init(void* _p, ae_state *_state)
{
    armijostate *p = (armijostate*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_init(&p->x, 0, DT_REAL, _state);
    ae_vector_init(&p->xbase, 0, DT_REAL, _state);
    ae_vector_init(&p->s, 0, DT_REAL, _state);
    _rcommstate_init(&p->rstate, _state);
}


void _armijostate_init_copy(void* _dst, void* _src, ae_state *_state)
{
    armijostate *dst = (armijostate*)_dst;
    armijostate *src = (armijostate*)_src;
    dst->needf = src->needf;
    ae_vector_init_copy(&dst->x, &src->x, _state);
    dst->f = src->f;
    dst->n = src->n;
    ae_vector_init_copy(&dst->xbase, &src->xbase, _state);
    ae_vector_init_copy(&dst->s, &src->s, _state);
    dst->stplen = src->stplen;
    dst->fcur = src->fcur;
    dst->stpmax = src->stpmax;
    dst->fmax = src->fmax;
    dst->nfev = src->nfev;
    dst->info = src->info;
    _rcommstate_init_copy(&dst->rstate, &src->rstate, _state);
}


void _armijostate_clear(void* _p)
{
    armijostate *p = (armijostate*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->xbase);
    ae_vector_clear(&p->s);
    _rcommstate_clear(&p->rstate);
}


void _armijostate_destroy(void* _p)
{
    armijostate *p = (armijostate*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_destroy(&p->x);
    ae_vector_destroy(&p->xbase);
    ae_vector_destroy(&p->s);
    _rcommstate_destroy(&p->rstate);
}




void findprimitiverootandinverse(ae_int_t n,
     ae_int_t* proot,
     ae_int_t* invproot,
     ae_state *_state)
{
    ae_int_t candroot;
    ae_int_t phin;
    ae_int_t q;
    ae_int_t f;
    ae_bool allnonone;
    ae_int_t x;
    ae_int_t lastx;
    ae_int_t y;
    ae_int_t lasty;
    ae_int_t a;
    ae_int_t b;
    ae_int_t t;
    ae_int_t n2;

    *proot = 0;
    *invproot = 0;

    ae_assert(n>=3, "FindPrimitiveRootAndInverse: N<3", _state);
    *proot = 0;
    *invproot = 0;
    
    /*
     * check that N is prime
     */
    ae_assert(ntheory_isprime(n, _state), "FindPrimitiveRoot: N is not prime", _state);
    
    /*
     * Because N is prime, Euler totient function is equal to N-1
     */
    phin = n-1;
    
    /*
     * Test different values of PRoot - from 2 to N-1.
     * One of these values MUST be primitive root.
     *
     * For testing we use algorithm from Wiki (Primitive root modulo n):
     * * compute phi(N)
     * * determine the different prime factors of phi(N), say p1, ..., pk
     * * for every element m of Zn*, compute m^(phi(N)/pi) mod N for i=1..k
     *   using a fast algorithm for modular exponentiation.
     * * a number m for which these k results are all different from 1 is a
     *   primitive root.
     */
    for(candroot=2; candroot<=n-1; candroot++)
    {
        
        /*
         * We have current candidate root in CandRoot.
         *
         * Scan different prime factors of PhiN. Here:
         * * F is a current candidate factor
         * * Q is a current quotient - amount which was left after dividing PhiN
         *   by all previous factors
         *
         * For each factor, perform test mentioned above.
         */
        q = phin;
        f = 2;
        allnonone = ae_true;
        while(q>1)
        {
            if( q%f==0 )
            {
                t = ntheory_modexp(candroot, phin/f, n, _state);
                if( t==1 )
                {
                    allnonone = ae_false;
                    break;
                }
                while(q%f==0)
                {
                    q = q/f;
                }
            }
            f = f+1;
        }
        if( allnonone )
        {
            *proot = candroot;
            break;
        }
    }
    ae_assert(*proot>=2, "FindPrimitiveRoot: internal error (root not found)", _state);
    
    /*
     * Use extended Euclidean algorithm to find multiplicative inverse of primitive root
     */
    x = 0;
    lastx = 1;
    y = 1;
    lasty = 0;
    a = *proot;
    b = n;
    while(b!=0)
    {
        q = a/b;
        t = a%b;
        a = b;
        b = t;
        t = lastx-q*x;
        lastx = x;
        x = t;
        t = lasty-q*y;
        lasty = y;
        y = t;
    }
    while(lastx<0)
    {
        lastx = lastx+n;
    }
    *invproot = lastx;
    
    /*
     * Check that it is safe to perform multiplication modulo N.
     * Check results for consistency.
     */
    n2 = (n-1)*(n-1);
    ae_assert(n2/(n-1)==n-1, "FindPrimitiveRoot: internal error", _state);
    ae_assert(*proot*(*invproot)/(*proot)==(*invproot), "FindPrimitiveRoot: internal error", _state);
    ae_assert(*proot*(*invproot)/(*invproot)==(*proot), "FindPrimitiveRoot: internal error", _state);
    ae_assert(*proot*(*invproot)%n==1, "FindPrimitiveRoot: internal error", _state);
}


static ae_bool ntheory_isprime(ae_int_t n, ae_state *_state)
{
    ae_int_t p;
    ae_bool result;


    result = ae_false;
    p = 2;
    while(p*p<=n)
    {
        if( n%p==0 )
        {
            return result;
        }
        p = p+1;
    }
    result = ae_true;
    return result;
}


static ae_int_t ntheory_modmul(ae_int_t a,
     ae_int_t b,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t t;
    double ra;
    double rb;
    ae_int_t result;


    ae_assert(a>=0&&a<n, "ModMul: A<0 or A>=N", _state);
    ae_assert(b>=0&&b<n, "ModMul: B<0 or B>=N", _state);
    
    /*
     * Base cases
     */
    ra = (double)(a);
    rb = (double)(b);
    if( b==0||a==0 )
    {
        result = 0;
        return result;
    }
    if( b==1||a==1 )
    {
        result = a*b;
        return result;
    }
    if( ae_fp_eq(ra*rb,(double)(a*b)) )
    {
        result = a*b%n;
        return result;
    }
    
    /*
     * Non-base cases
     */
    if( b%2==0 )
    {
        
        /*
         * A*B = (A*(B/2)) * 2
         *
         * Product T=A*(B/2) is calculated recursively, product T*2 is
         * calculated as follows:
         * * result:=T-N
         * * result:=result+T
         * * if result<0 then result:=result+N
         *
         * In case integer result overflows, we generate exception
         */
        t = ntheory_modmul(a, b/2, n, _state);
        result = t-n;
        result = result+t;
        if( result<0 )
        {
            result = result+n;
        }
    }
    else
    {
        
        /*
         * A*B = (A*(B div 2)) * 2 + A
         *
         * Product T=A*(B/2) is calculated recursively, product T*2 is
         * calculated as follows:
         * * result:=T-N
         * * result:=result+T
         * * if result<0 then result:=result+N
         *
         * In case integer result overflows, we generate exception
         */
        t = ntheory_modmul(a, b/2, n, _state);
        result = t-n;
        result = result+t;
        if( result<0 )
        {
            result = result+n;
        }
        result = result-n;
        result = result+a;
        if( result<0 )
        {
            result = result+n;
        }
    }
    return result;
}


static ae_int_t ntheory_modexp(ae_int_t a,
     ae_int_t b,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t t;
    ae_int_t result;


    ae_assert(a>=0&&a<n, "ModExp: A<0 or A>=N", _state);
    ae_assert(b>=0, "ModExp: B<0", _state);
    
    /*
     * Base cases
     */
    if( b==0 )
    {
        result = 1;
        return result;
    }
    if( b==1 )
    {
        result = a;
        return result;
    }
    
    /*
     * Non-base cases
     */
    if( b%2==0 )
    {
        t = ntheory_modmul(a, a, n, _state);
        result = ntheory_modexp(t, b/2, n, _state);
    }
    else
    {
        t = ntheory_modmul(a, a, n, _state);
        result = ntheory_modexp(t, b/2, n, _state);
        result = ntheory_modmul(result, a, n, _state);
    }
    return result;
}




/*************************************************************************
This subroutine generates FFT plan for K complex FFT's with length N each.

INPUT PARAMETERS:
    N           -   FFT length (in complex numbers), N>=1
    K           -   number of repetitions, K>=1
    
OUTPUT PARAMETERS:
    Plan        -   plan

  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
void ftcomplexfftplan(ae_int_t n,
     ae_int_t k,
     fasttransformplan* plan,
     ae_state *_state)
{
    ae_frame _frame_block;
    srealarray bluesteinbuf;
    ae_int_t rowptr;
    ae_int_t bluesteinsize;
    ae_int_t precrptr;
    ae_int_t preciptr;
    ae_int_t precrsize;
    ae_int_t precisize;

    ae_frame_make(_state, &_frame_block);
    _fasttransformplan_clear(plan);
    _srealarray_init(&bluesteinbuf, _state);

    
    /*
     * Initial check for parameters
     */
    ae_assert(n>0, "FTComplexFFTPlan: N<=0", _state);
    ae_assert(k>0, "FTComplexFFTPlan: K<=0", _state);
    
    /*
     * Determine required sizes of precomputed real and integer
     * buffers. This stage of code is highly dependent on internals
     * of FTComplexFFTPlanRec() and must be kept synchronized with
     * possible changes in internals of plan generation function.
     *
     * Buffer size is determined as follows:
     * * N is factorized
     * * we factor out anything which is less or equal to MaxRadix
     * * prime factor F>RaderThreshold requires 4*FTBaseFindSmooth(2*F-1)
     *   real entries to store precomputed Quantities for Bluestein's
     *   transformation
     * * prime factor F<=RaderThreshold does NOT require
     *   precomputed storage
     */
    precrsize = 0;
    precisize = 0;
    ftbase_ftdeterminespacerequirements(n, &precrsize, &precisize, _state);
    if( precrsize>0 )
    {
        ae_vector_set_length(&plan->precr, precrsize, _state);
    }
    if( precisize>0 )
    {
        ae_vector_set_length(&plan->preci, precisize, _state);
    }
    
    /*
     * Generate plan
     */
    rowptr = 0;
    precrptr = 0;
    preciptr = 0;
    bluesteinsize = 1;
    ae_vector_set_length(&plan->buffer, 2*n*k, _state);
    ftbase_ftcomplexfftplanrec(n, k, ae_true, ae_true, &rowptr, &bluesteinsize, &precrptr, &preciptr, plan, _state);
    ae_vector_set_length(&bluesteinbuf.val, bluesteinsize, _state);
    ae_shared_pool_set_seed(&plan->bluesteinpool, &bluesteinbuf, sizeof(bluesteinbuf), _srealarray_init, _srealarray_init_copy, _srealarray_destroy, _state);
    
    /*
     * Check that actual amount of precomputed space used by transformation
     * plan is EXACTLY equal to amount of space allocated by us.
     */
    ae_assert(precrptr==precrsize, "FTComplexFFTPlan: internal error (PrecRPtr<>PrecRSize)", _state);
    ae_assert(preciptr==precisize, "FTComplexFFTPlan: internal error (PrecRPtr<>PrecRSize)", _state);
    ae_frame_leave(_state);
}


/*************************************************************************
This subroutine applies transformation plan to input/output array A.

INPUT PARAMETERS:
    Plan        -   transformation plan
    A           -   array, must be large enough for plan to work
    OffsA       -   offset of the subarray to process
    RepCnt      -   repetition count (transformation is repeatedly applied
                    to subsequent subarrays)
    
OUTPUT PARAMETERS:
    Plan        -   plan (temporary buffers can be modified, plan itself
                    is unchanged and can be reused)
    A           -   transformed array

  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
void ftapplyplan(fasttransformplan* plan,
     /* Real    */ ae_vector* a,
     ae_int_t offsa,
     ae_int_t repcnt,
     ae_state *_state)
{
    ae_int_t plansize;
    ae_int_t i;


    plansize = plan->entries.ptr.pp_int[0][ftbase_coloperandscnt]*plan->entries.ptr.pp_int[0][ftbase_coloperandsize]*plan->entries.ptr.pp_int[0][ftbase_colmicrovectorsize];
    for(i=0; i<=repcnt-1; i++)
    {
        ftbase_ftapplysubplan(plan, 0, a, offsa+plansize*i, 0, &plan->buffer, 1, _state);
    }
}


/*************************************************************************
Returns good factorization N=N1*N2.

Usually N1<=N2 (but not always - small N's may be exception).
if N1<>1 then N2<>1.

Factorization is chosen depending on task type and codelets we have.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbasefactorize(ae_int_t n,
     ae_int_t tasktype,
     ae_int_t* n1,
     ae_int_t* n2,
     ae_state *_state)
{
    ae_int_t j;

    *n1 = 0;
    *n2 = 0;

    *n1 = 0;
    *n2 = 0;
    
    /*
     * try to find good codelet
     */
    if( *n1*(*n2)!=n )
    {
        for(j=ftbase_ftbasecodeletrecommended; j>=2; j--)
        {
            if( n%j==0 )
            {
                *n1 = j;
                *n2 = n/j;
                break;
            }
        }
    }
    
    /*
     * try to factorize N
     */
    if( *n1*(*n2)!=n )
    {
        for(j=ftbase_ftbasecodeletrecommended+1; j<=n-1; j++)
        {
            if( n%j==0 )
            {
                *n1 = j;
                *n2 = n/j;
                break;
            }
        }
    }
    
    /*
     * looks like N is prime :(
     */
    if( *n1*(*n2)!=n )
    {
        *n1 = 1;
        *n2 = n;
    }
    
    /*
     * normalize
     */
    if( *n2==1&&*n1!=1 )
    {
        *n2 = *n1;
        *n1 = 1;
    }
}


/*************************************************************************
Is number smooth?

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool ftbaseissmooth(ae_int_t n, ae_state *_state)
{
    ae_int_t i;
    ae_bool result;


    for(i=2; i<=ftbase_ftbasemaxsmoothfactor; i++)
    {
        while(n%i==0)
        {
            n = n/i;
        }
    }
    result = n==1;
    return result;
}


/*************************************************************************
Returns smallest smooth (divisible only by 2, 3, 5) number that is greater
than or equal to max(N,2)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
ae_int_t ftbasefindsmooth(ae_int_t n, ae_state *_state)
{
    ae_int_t best;
    ae_int_t result;


    best = 2;
    while(best<n)
    {
        best = 2*best;
    }
    ftbase_ftbasefindsmoothrec(n, 1, 2, &best, _state);
    result = best;
    return result;
}


/*************************************************************************
Returns  smallest  smooth  (divisible only by 2, 3, 5) even number that is
greater than or equal to max(N,2)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
ae_int_t ftbasefindsmootheven(ae_int_t n, ae_state *_state)
{
    ae_int_t best;
    ae_int_t result;


    best = 2;
    while(best<n)
    {
        best = 2*best;
    }
    ftbase_ftbasefindsmoothrec(n, 2, 2, &best, _state);
    result = best;
    return result;
}


/*************************************************************************
Returns estimate of FLOP count for the FFT.

It is only an estimate based on operations count for the PERFECT FFT
and relative inefficiency of the algorithm actually used.

N should be power of 2, estimates are badly wrong for non-power-of-2 N's.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
double ftbasegetflopestimate(ae_int_t n, ae_state *_state)
{
    double result;


    result = ftbase_ftbaseinefficiencyfactor*(4*n*ae_log((double)(n), _state)/ae_log((double)(2), _state)-6*n+8);
    return result;
}


/*************************************************************************
This function returns EXACT estimate of the space requirements for N-point
FFT. Internals of this function are highly dependent on details of different
FFTs employed by this unit, so every time algorithm is changed this function
has to be rewritten.

INPUT PARAMETERS:
    N           -   transform length
    PrecRSize   -   must be set to zero
    PrecISize   -   must be set to zero
    
OUTPUT PARAMETERS:
    PrecRSize   -   number of real temporaries required for transformation
    PrecISize   -   number of integer temporaries required for transformation

    
  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftdeterminespacerequirements(ae_int_t n,
     ae_int_t* precrsize,
     ae_int_t* precisize,
     ae_state *_state)
{
    ae_int_t ncur;
    ae_int_t f;
    ae_int_t i;


    
    /*
     * Determine required sizes of precomputed real and integer
     * buffers. This stage of code is highly dependent on internals
     * of FTComplexFFTPlanRec() and must be kept synchronized with
     * possible changes in internals of plan generation function.
     *
     * Buffer size is determined as follows:
     * * N is factorized
     * * we factor out anything which is less or equal to MaxRadix
     * * prime factor F>RaderThreshold requires 4*FTBaseFindSmooth(2*F-1)
     *   real entries to store precomputed Quantities for Bluestein's
     *   transformation
     * * prime factor F<=RaderThreshold requires 2*(F-1)+ESTIMATE(F-1)
     *   precomputed storage
     */
    ncur = n;
    for(i=2; i<=ftbase_maxradix; i++)
    {
        while(ncur%i==0)
        {
            ncur = ncur/i;
        }
    }
    f = 2;
    while(f<=ncur)
    {
        while(ncur%f==0)
        {
            if( f>ftbase_raderthreshold )
            {
                *precrsize = *precrsize+4*ftbasefindsmooth(2*f-1, _state);
            }
            else
            {
                *precrsize = *precrsize+2*(f-1);
                ftbase_ftdeterminespacerequirements(f-1, precrsize, precisize, _state);
            }
            ncur = ncur/f;
        }
        f = f+1;
    }
}


/*************************************************************************
Recurrent function called by FTComplexFFTPlan() and other functions. It
recursively builds transformation plan

INPUT PARAMETERS:
    N           -   FFT length (in complex numbers), N>=1
    K           -   number of repetitions, K>=1
    ChildPlan   -   if True, plan generator inserts OpStart/opEnd in the
                    plan header/footer.
    TopmostPlan -   if True, plan generator assumes that it is topmost plan:
                    * it may use global buffer for transpositions
                    and there is no other plan which executes in parallel
    RowPtr      -   index which points to past-the-last entry generated so far
    BluesteinSize-  amount of storage (in real numbers) required for Bluestein buffer
    PrecRPtr    -   pointer to unused part of precomputed real buffer (Plan.PrecR):
                    * when this function stores some data to precomputed buffer,
                      it advances pointer.
                    * it is responsibility of the function to assert that
                      Plan.PrecR has enough space to store data before actually
                      writing to buffer.
                    * it is responsibility of the caller to allocate enough
                      space before calling this function
    PrecIPtr    -   pointer to unused part of precomputed integer buffer (Plan.PrecI):
                    * when this function stores some data to precomputed buffer,
                      it advances pointer.
                    * it is responsibility of the function to assert that
                      Plan.PrecR has enough space to store data before actually
                      writing to buffer.
                    * it is responsibility of the caller to allocate enough
                      space before calling this function
    Plan        -   plan (generated so far)
    
OUTPUT PARAMETERS:
    RowPtr      -   updated pointer (advanced by number of entries generated
                    by function)
    BluesteinSize-  updated amount
                    (may be increased, but may never be decreased)
        
NOTE: in case TopmostPlan is True, ChildPlan is also must be True.
    
  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftcomplexfftplanrec(ae_int_t n,
     ae_int_t k,
     ae_bool childplan,
     ae_bool topmostplan,
     ae_int_t* rowptr,
     ae_int_t* bluesteinsize,
     ae_int_t* precrptr,
     ae_int_t* preciptr,
     fasttransformplan* plan,
     ae_state *_state)
{
    ae_frame _frame_block;
    srealarray localbuf;
    ae_int_t m;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t gq;
    ae_int_t giq;
    ae_int_t row0;
    ae_int_t row1;
    ae_int_t row2;
    ae_int_t row3;

    ae_frame_make(_state, &_frame_block);
    _srealarray_init(&localbuf, _state);

    ae_assert(n>0, "FTComplexFFTPlan: N<=0", _state);
    ae_assert(k>0, "FTComplexFFTPlan: K<=0", _state);
    ae_assert(!topmostplan||childplan, "FTComplexFFTPlan: ChildPlan is inconsistent with TopmostPlan", _state);
    
    /*
     * Try to generate "topmost" plan
     */
    if( topmostplan&&n>ftbase_recursivethreshold )
    {
        ftbase_ftfactorize(n, ae_false, &n1, &n2, _state);
        if( n1*n2==0 )
        {
            
            /*
             * Handle prime-factor FFT with Bluestein's FFT.
             * Determine size of Bluestein's buffer.
             */
            m = ftbasefindsmooth(2*n-1, _state);
            *bluesteinsize = ae_maxint(2*m, *bluesteinsize, _state);
            
            /*
             * Generate plan
             */
            ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
            ftbase_ftpushentry4(plan, rowptr, ftbase_opbluesteinsfft, k, n, 2, m, 2, *precrptr, 0, _state);
            row0 = *rowptr;
            ftbase_ftpushentry(plan, rowptr, ftbase_opjmp, 0, 0, 0, 0, _state);
            ftbase_ftcomplexfftplanrec(m, 1, ae_true, ae_true, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
            row1 = *rowptr;
            plan->entries.ptr.pp_int[row0][ftbase_colparam0] = row1-row0;
            ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
            
            /*
             * Fill precomputed buffer
             */
            ftbase_ftprecomputebluesteinsfft(n, m, &plan->precr, *precrptr, _state);
            
            /*
             * Update pointer to the precomputed area
             */
            *precrptr = *precrptr+4*m;
        }
        else
        {
            
            /*
             * Handle composite FFT with recursive Cooley-Tukey which
             * uses global buffer instead of local one.
             */
            ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
            ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
            row0 = *rowptr;
            ftbase_ftpushentry2(plan, rowptr, ftbase_opparallelcall, k*n2, n1, 2, 0, ftbase_ftoptimisticestimate(n, _state), _state);
            ftbase_ftpushentry(plan, rowptr, ftbase_opcomplexfftfactors, k, n, 2, n1, _state);
            ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n2, _state);
            row2 = *rowptr;
            ftbase_ftpushentry2(plan, rowptr, ftbase_opparallelcall, k*n1, n2, 2, 0, ftbase_ftoptimisticestimate(n, _state), _state);
            ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
            ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
            row1 = *rowptr;
            ftbase_ftcomplexfftplanrec(n1, 1, ae_true, ae_false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
            plan->entries.ptr.pp_int[row0][ftbase_colparam0] = row1-row0;
            row3 = *rowptr;
            ftbase_ftcomplexfftplanrec(n2, 1, ae_true, ae_false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
            plan->entries.ptr.pp_int[row2][ftbase_colparam0] = row3-row2;
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Prepare "non-topmost" plan:
     * * calculate factorization
     * * use local (shared) buffer
     * * update buffer size - ANY plan will need at least
     *   2*N temporaries, additional requirements can be
     *   applied later
     */
    ftbase_ftfactorize(n, ae_false, &n1, &n2, _state);
    
    /*
     * Handle FFT's with N1*N2=0: either small-N or prime-factor
     */
    if( n1*n2==0 )
    {
        if( n<=ftbase_maxradix )
        {
            
            /*
             * Small-N FFT
             */
            if( childplan )
            {
                ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
            }
            ftbase_ftpushentry(plan, rowptr, ftbase_opcomplexcodeletfft, k, n, 2, 0, _state);
            if( childplan )
            {
                ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
            }
            ae_frame_leave(_state);
            return;
        }
        if( n<=ftbase_raderthreshold )
        {
            
            /*
             * Handle prime-factor FFT's with Rader's FFT
             */
            m = n-1;
            if( childplan )
            {
                ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
            }
            findprimitiverootandinverse(n, &gq, &giq, _state);
            ftbase_ftpushentry4(plan, rowptr, ftbase_opradersfft, k, n, 2, 2, gq, giq, *precrptr, _state);
            ftbase_ftprecomputeradersfft(n, gq, giq, &plan->precr, *precrptr, _state);
            *precrptr = *precrptr+2*(n-1);
            row0 = *rowptr;
            ftbase_ftpushentry(plan, rowptr, ftbase_opjmp, 0, 0, 0, 0, _state);
            ftbase_ftcomplexfftplanrec(m, 1, ae_true, ae_false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
            row1 = *rowptr;
            plan->entries.ptr.pp_int[row0][ftbase_colparam0] = row1-row0;
            if( childplan )
            {
                ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
            }
        }
        else
        {
            
            /*
             * Handle prime-factor FFT's with Bluestein's FFT
             */
            m = ftbasefindsmooth(2*n-1, _state);
            *bluesteinsize = ae_maxint(2*m, *bluesteinsize, _state);
            if( childplan )
            {
                ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
            }
            ftbase_ftpushentry4(plan, rowptr, ftbase_opbluesteinsfft, k, n, 2, m, 2, *precrptr, 0, _state);
            ftbase_ftprecomputebluesteinsfft(n, m, &plan->precr, *precrptr, _state);
            *precrptr = *precrptr+4*m;
            row0 = *rowptr;
            ftbase_ftpushentry(plan, rowptr, ftbase_opjmp, 0, 0, 0, 0, _state);
            ftbase_ftcomplexfftplanrec(m, 1, ae_true, ae_false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
            row1 = *rowptr;
            plan->entries.ptr.pp_int[row0][ftbase_colparam0] = row1-row0;
            if( childplan )
            {
                ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
            }
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Handle Cooley-Tukey FFT with small N1
     */
    if( n1<=ftbase_maxradix )
    {
        
        /*
         * Specialized transformation for small N1:
         * * N2 short inplace FFT's, each N1-point, with integrated twiddle factors
         * * N1 long FFT's
         * * final transposition
         */
        if( childplan )
        {
            ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
        }
        ftbase_ftpushentry(plan, rowptr, ftbase_opcomplexcodelettwfft, k, n1, 2*n2, 0, _state);
        ftbase_ftcomplexfftplanrec(n2, k*n1, ae_false, ae_false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
        ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
        if( childplan )
        {
            ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Handle general Cooley-Tukey FFT, either "flat" or "recursive"
     */
    if( n<=ftbase_recursivethreshold )
    {
        
        /*
         * General code for large N1/N2, "flat" version without explicit recurrence
         * (nested subplans are inserted directly into the body of the plan)
         */
        if( childplan )
        {
            ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
        }
        ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
        ftbase_ftcomplexfftplanrec(n1, k*n2, ae_false, ae_false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
        ftbase_ftpushentry(plan, rowptr, ftbase_opcomplexfftfactors, k, n, 2, n1, _state);
        ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n2, _state);
        ftbase_ftcomplexfftplanrec(n2, k*n1, ae_false, ae_false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
        ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
        if( childplan )
        {
            ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
        }
    }
    else
    {
        
        /*
         * General code for large N1/N2, "recursive" version - nested subplans
         * are separated from the plan body.
         *
         * Generate parent plan.
         */
        if( childplan )
        {
            ftbase_ftpushentry2(plan, rowptr, ftbase_opstart, k, n, 2, -1, ftbase_ftoptimisticestimate(n, _state), _state);
        }
        ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
        row0 = *rowptr;
        ftbase_ftpushentry2(plan, rowptr, ftbase_opparallelcall, k*n2, n1, 2, 0, ftbase_ftoptimisticestimate(n, _state), _state);
        ftbase_ftpushentry(plan, rowptr, ftbase_opcomplexfftfactors, k, n, 2, n1, _state);
        ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n2, _state);
        row2 = *rowptr;
        ftbase_ftpushentry2(plan, rowptr, ftbase_opparallelcall, k*n1, n2, 2, 0, ftbase_ftoptimisticestimate(n, _state), _state);
        ftbase_ftpushentry(plan, rowptr, ftbase_opcomplextranspose, k, n, 2, n1, _state);
        if( childplan )
        {
            ftbase_ftpushentry(plan, rowptr, ftbase_opend, k, n, 2, 0, _state);
        }
        
        /*
         * Generate child subplans, insert refence to parent plans
         */
        row1 = *rowptr;
        ftbase_ftcomplexfftplanrec(n1, 1, ae_true, ae_false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
        plan->entries.ptr.pp_int[row0][ftbase_colparam0] = row1-row0;
        row3 = *rowptr;
        ftbase_ftcomplexfftplanrec(n2, 1, ae_true, ae_false, rowptr, bluesteinsize, precrptr, preciptr, plan, _state);
        plan->entries.ptr.pp_int[row2][ftbase_colparam0] = row3-row2;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This function pushes one more entry to the plan. It resizes Entries matrix
if needed.

INPUT PARAMETERS:
    Plan        -   plan (generated so far)
    RowPtr      -   index which points to past-the-last entry generated so far
    EType       -   entry type
    EOpCnt      -   operands count
    EOpSize     -   operand size
    EMcvSize    -   microvector size
    EParam0     -   parameter 0
    
OUTPUT PARAMETERS:
    Plan        -   updated plan    
    RowPtr      -   updated pointer

NOTE: Param1 is set to -1.
    
  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftpushentry(fasttransformplan* plan,
     ae_int_t* rowptr,
     ae_int_t etype,
     ae_int_t eopcnt,
     ae_int_t eopsize,
     ae_int_t emcvsize,
     ae_int_t eparam0,
     ae_state *_state)
{


    ftbase_ftpushentry2(plan, rowptr, etype, eopcnt, eopsize, emcvsize, eparam0, -1, _state);
}


/*************************************************************************
Same as FTPushEntry(), but sets Param0 AND Param1.
This function pushes one more entry to the plan. It resized Entries matrix
if needed.

INPUT PARAMETERS:
    Plan        -   plan (generated so far)
    RowPtr      -   index which points to past-the-last entry generated so far
    EType       -   entry type
    EOpCnt      -   operands count
    EOpSize     -   operand size
    EMcvSize    -   microvector size
    EParam0     -   parameter 0
    EParam1     -   parameter 1
    
OUTPUT PARAMETERS:
    Plan        -   updated plan    
    RowPtr      -   updated pointer

  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftpushentry2(fasttransformplan* plan,
     ae_int_t* rowptr,
     ae_int_t etype,
     ae_int_t eopcnt,
     ae_int_t eopsize,
     ae_int_t emcvsize,
     ae_int_t eparam0,
     ae_int_t eparam1,
     ae_state *_state)
{


    if( *rowptr>=plan->entries.rows )
    {
        imatrixresize(&plan->entries, ae_maxint(2*plan->entries.rows, 1, _state), ftbase_colscnt, _state);
    }
    plan->entries.ptr.pp_int[*rowptr][ftbase_coltype] = etype;
    plan->entries.ptr.pp_int[*rowptr][ftbase_coloperandscnt] = eopcnt;
    plan->entries.ptr.pp_int[*rowptr][ftbase_coloperandsize] = eopsize;
    plan->entries.ptr.pp_int[*rowptr][ftbase_colmicrovectorsize] = emcvsize;
    plan->entries.ptr.pp_int[*rowptr][ftbase_colparam0] = eparam0;
    plan->entries.ptr.pp_int[*rowptr][ftbase_colparam1] = eparam1;
    plan->entries.ptr.pp_int[*rowptr][ftbase_colparam2] = 0;
    plan->entries.ptr.pp_int[*rowptr][ftbase_colparam3] = 0;
    *rowptr = *rowptr+1;
}


/*************************************************************************
Same as FTPushEntry(), but sets Param0, Param1, Param2 and Param3.
This function pushes one more entry to the plan. It resized Entries matrix
if needed.

INPUT PARAMETERS:
    Plan        -   plan (generated so far)
    RowPtr      -   index which points to past-the-last entry generated so far
    EType       -   entry type
    EOpCnt      -   operands count
    EOpSize     -   operand size
    EMcvSize    -   microvector size
    EParam0     -   parameter 0
    EParam1     -   parameter 1
    EParam2     -   parameter 2
    EParam3     -   parameter 3
    
OUTPUT PARAMETERS:
    Plan        -   updated plan    
    RowPtr      -   updated pointer

  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftpushentry4(fasttransformplan* plan,
     ae_int_t* rowptr,
     ae_int_t etype,
     ae_int_t eopcnt,
     ae_int_t eopsize,
     ae_int_t emcvsize,
     ae_int_t eparam0,
     ae_int_t eparam1,
     ae_int_t eparam2,
     ae_int_t eparam3,
     ae_state *_state)
{


    if( *rowptr>=plan->entries.rows )
    {
        imatrixresize(&plan->entries, ae_maxint(2*plan->entries.rows, 1, _state), ftbase_colscnt, _state);
    }
    plan->entries.ptr.pp_int[*rowptr][ftbase_coltype] = etype;
    plan->entries.ptr.pp_int[*rowptr][ftbase_coloperandscnt] = eopcnt;
    plan->entries.ptr.pp_int[*rowptr][ftbase_coloperandsize] = eopsize;
    plan->entries.ptr.pp_int[*rowptr][ftbase_colmicrovectorsize] = emcvsize;
    plan->entries.ptr.pp_int[*rowptr][ftbase_colparam0] = eparam0;
    plan->entries.ptr.pp_int[*rowptr][ftbase_colparam1] = eparam1;
    plan->entries.ptr.pp_int[*rowptr][ftbase_colparam2] = eparam2;
    plan->entries.ptr.pp_int[*rowptr][ftbase_colparam3] = eparam3;
    *rowptr = *rowptr+1;
}


/*************************************************************************
This subroutine applies subplan to input/output array A.

INPUT PARAMETERS:
    Plan        -   transformation plan
    SubPlan     -   subplan index
    A           -   array, must be large enough for plan to work
    ABase       -   base offset in array A, this value points to start of
                    subarray whose length is equal to length of the plan
    AOffset     -   offset with respect to ABase, 0<=AOffset<PlanLength.
                    This is an offset within large PlanLength-subarray of
                    the chunk to process.
    Buf         -   temporary buffer whose length is equal to plan length
                    (without taking into account RepCnt) or larger.
    OffsBuf     -   offset in the buffer array
    RepCnt      -   repetition count (transformation is repeatedly applied
                    to subsequent subarrays)
    
OUTPUT PARAMETERS:
    Plan        -   plan (temporary buffers can be modified, plan itself
                    is unchanged and can be reused)
    A           -   transformed array

  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftapplysubplan(fasttransformplan* plan,
     ae_int_t subplan,
     /* Real    */ ae_vector* a,
     ae_int_t abase,
     ae_int_t aoffset,
     /* Real    */ ae_vector* buf,
     ae_int_t repcnt,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t rowidx;
    ae_int_t i;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t operation;
    ae_int_t operandscnt;
    ae_int_t operandsize;
    ae_int_t microvectorsize;
    ae_int_t param0;
    ae_int_t param1;
    ae_int_t parentsize;
    ae_int_t childsize;
    ae_int_t chunksize;
    ae_int_t lastchunksize;
    srealarray *bufa;
    ae_smart_ptr _bufa;
    srealarray *bufb;
    ae_smart_ptr _bufb;
    srealarray *bufc;
    ae_smart_ptr _bufc;
    srealarray *bufd;
    ae_smart_ptr _bufd;

    ae_frame_make(_state, &_frame_block);
    ae_smart_ptr_init(&_bufa, (void**)&bufa, _state);
    ae_smart_ptr_init(&_bufb, (void**)&bufb, _state);
    ae_smart_ptr_init(&_bufc, (void**)&bufc, _state);
    ae_smart_ptr_init(&_bufd, (void**)&bufd, _state);

    ae_assert(plan->entries.ptr.pp_int[subplan][ftbase_coltype]==ftbase_opstart, "FTApplySubPlan: incorrect subplan header", _state);
    rowidx = subplan+1;
    while(plan->entries.ptr.pp_int[rowidx][ftbase_coltype]!=ftbase_opend)
    {
        operation = plan->entries.ptr.pp_int[rowidx][ftbase_coltype];
        operandscnt = repcnt*plan->entries.ptr.pp_int[rowidx][ftbase_coloperandscnt];
        operandsize = plan->entries.ptr.pp_int[rowidx][ftbase_coloperandsize];
        microvectorsize = plan->entries.ptr.pp_int[rowidx][ftbase_colmicrovectorsize];
        param0 = plan->entries.ptr.pp_int[rowidx][ftbase_colparam0];
        param1 = plan->entries.ptr.pp_int[rowidx][ftbase_colparam1];
        touchint(&param1, _state);
        
        /*
         * Process "jump" operation
         */
        if( operation==ftbase_opjmp )
        {
            rowidx = rowidx+plan->entries.ptr.pp_int[rowidx][ftbase_colparam0];
            continue;
        }
        
        /*
         * Process "parallel call" operation:
         * * we perform initial check for consistency between parent and child plans
         * * we call FTSplitAndApplyParallelPlan(), which splits parallel plan into
         *   several parallel tasks
         */
        if( operation==ftbase_opparallelcall )
        {
            parentsize = operandsize*microvectorsize;
            childsize = plan->entries.ptr.pp_int[rowidx+param0][ftbase_coloperandscnt]*plan->entries.ptr.pp_int[rowidx+param0][ftbase_coloperandsize]*plan->entries.ptr.pp_int[rowidx+param0][ftbase_colmicrovectorsize];
            ae_assert(plan->entries.ptr.pp_int[rowidx+param0][ftbase_coltype]==ftbase_opstart, "FTApplySubPlan: incorrect child subplan header", _state);
            ae_assert(parentsize==childsize, "FTApplySubPlan: incorrect child subplan header", _state);
            chunksize = ae_maxint(ftbase_recursivethreshold/childsize, 1, _state);
            lastchunksize = operandscnt%chunksize;
            if( lastchunksize==0 )
            {
                lastchunksize = chunksize;
            }
            i = 0;
            while(i<operandscnt)
            {
                chunksize = ae_minint(chunksize, operandscnt-i, _state);
                ftbase_ftapplysubplan(plan, rowidx+param0, a, abase, aoffset+i*childsize, buf, chunksize, _state);
                i = i+chunksize;
            }
            rowidx = rowidx+1;
            continue;
        }
        
        /*
         * Process "reference complex FFT" operation
         */
        if( operation==ftbase_opcomplexreffft )
        {
            ftbase_ftapplycomplexreffft(a, abase+aoffset, operandscnt, operandsize, microvectorsize, buf, _state);
            rowidx = rowidx+1;
            continue;
        }
        
        /*
         * Process "codelet FFT" operation
         */
        if( operation==ftbase_opcomplexcodeletfft )
        {
            ftbase_ftapplycomplexcodeletfft(a, abase+aoffset, operandscnt, operandsize, microvectorsize, _state);
            rowidx = rowidx+1;
            continue;
        }
        
        /*
         * Process "integrated codelet FFT" operation
         */
        if( operation==ftbase_opcomplexcodelettwfft )
        {
            ftbase_ftapplycomplexcodelettwfft(a, abase+aoffset, operandscnt, operandsize, microvectorsize, _state);
            rowidx = rowidx+1;
            continue;
        }
        
        /*
         * Process Bluestein's FFT operation
         */
        if( operation==ftbase_opbluesteinsfft )
        {
            ae_assert(microvectorsize==2, "FTApplySubPlan: microvectorsize!=2 for Bluesteins FFT", _state);
            ae_shared_pool_retrieve(&plan->bluesteinpool, &_bufa, _state);
            ae_shared_pool_retrieve(&plan->bluesteinpool, &_bufb, _state);
            ae_shared_pool_retrieve(&plan->bluesteinpool, &_bufc, _state);
            ae_shared_pool_retrieve(&plan->bluesteinpool, &_bufd, _state);
            ftbase_ftbluesteinsfft(plan, a, abase, aoffset, operandscnt, operandsize, plan->entries.ptr.pp_int[rowidx][ftbase_colparam0], plan->entries.ptr.pp_int[rowidx][ftbase_colparam2], rowidx+plan->entries.ptr.pp_int[rowidx][ftbase_colparam1], &bufa->val, &bufb->val, &bufc->val, &bufd->val, _state);
            ae_shared_pool_recycle(&plan->bluesteinpool, &_bufa, _state);
            ae_shared_pool_recycle(&plan->bluesteinpool, &_bufb, _state);
            ae_shared_pool_recycle(&plan->bluesteinpool, &_bufc, _state);
            ae_shared_pool_recycle(&plan->bluesteinpool, &_bufd, _state);
            rowidx = rowidx+1;
            continue;
        }
        
        /*
         * Process Rader's FFT
         */
        if( operation==ftbase_opradersfft )
        {
            ftbase_ftradersfft(plan, a, abase, aoffset, operandscnt, operandsize, rowidx+plan->entries.ptr.pp_int[rowidx][ftbase_colparam0], plan->entries.ptr.pp_int[rowidx][ftbase_colparam1], plan->entries.ptr.pp_int[rowidx][ftbase_colparam2], plan->entries.ptr.pp_int[rowidx][ftbase_colparam3], buf, _state);
            rowidx = rowidx+1;
            continue;
        }
        
        /*
         * Process "complex twiddle factors" operation
         */
        if( operation==ftbase_opcomplexfftfactors )
        {
            ae_assert(microvectorsize==2, "FTApplySubPlan: MicrovectorSize<>1", _state);
            n1 = plan->entries.ptr.pp_int[rowidx][ftbase_colparam0];
            n2 = operandsize/n1;
            for(i=0; i<=operandscnt-1; i++)
            {
                ftbase_ffttwcalc(a, abase+aoffset+i*operandsize*2, n1, n2, _state);
            }
            rowidx = rowidx+1;
            continue;
        }
        
        /*
         * Process "complex transposition" operation
         */
        if( operation==ftbase_opcomplextranspose )
        {
            ae_assert(microvectorsize==2, "FTApplySubPlan: MicrovectorSize<>1", _state);
            n1 = plan->entries.ptr.pp_int[rowidx][ftbase_colparam0];
            n2 = operandsize/n1;
            for(i=0; i<=operandscnt-1; i++)
            {
                ftbase_internalcomplexlintranspose(a, n1, n2, abase+aoffset+i*operandsize*2, buf, _state);
            }
            rowidx = rowidx+1;
            continue;
        }
        
        /*
         * Error
         */
        ae_assert(ae_false, "FTApplySubPlan: unexpected plan type", _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This subroutine applies complex reference FFT to input/output array A.

VERY SLOW OPERATION, do not use it in real life plans :)

INPUT PARAMETERS:
    A           -   array, must be large enough for plan to work
    Offs        -   offset of the subarray to process
    OperandsCnt -   operands count (see description of FastTransformPlan)
    OperandSize -   operand size (see description of FastTransformPlan)
    MicrovectorSize-microvector size (see description of FastTransformPlan)
    Buf         -   temporary array, must be at least OperandsCnt*OperandSize*MicrovectorSize
    
OUTPUT PARAMETERS:
    A           -   transformed array

  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftapplycomplexreffft(/* Real    */ ae_vector* a,
     ae_int_t offs,
     ae_int_t operandscnt,
     ae_int_t operandsize,
     ae_int_t microvectorsize,
     /* Real    */ ae_vector* buf,
     ae_state *_state)
{
    ae_int_t opidx;
    ae_int_t i;
    ae_int_t k;
    double hre;
    double him;
    double c;
    double s;
    double re;
    double im;
    ae_int_t n;


    ae_assert(operandscnt>=1, "FTApplyComplexRefFFT: OperandsCnt<1", _state);
    ae_assert(operandsize>=1, "FTApplyComplexRefFFT: OperandSize<1", _state);
    ae_assert(microvectorsize==2, "FTApplyComplexRefFFT: MicrovectorSize<>2", _state);
    n = operandsize;
    for(opidx=0; opidx<=operandscnt-1; opidx++)
    {
        for(i=0; i<=n-1; i++)
        {
            hre = (double)(0);
            him = (double)(0);
            for(k=0; k<=n-1; k++)
            {
                re = a->ptr.p_double[offs+opidx*operandsize*2+2*k+0];
                im = a->ptr.p_double[offs+opidx*operandsize*2+2*k+1];
                c = ae_cos(-2*ae_pi*k*i/n, _state);
                s = ae_sin(-2*ae_pi*k*i/n, _state);
                hre = hre+c*re-s*im;
                him = him+c*im+s*re;
            }
            buf->ptr.p_double[2*i+0] = hre;
            buf->ptr.p_double[2*i+1] = him;
        }
        for(i=0; i<=operandsize*2-1; i++)
        {
            a->ptr.p_double[offs+opidx*operandsize*2+i] = buf->ptr.p_double[i];
        }
    }
}


/*************************************************************************
This subroutine applies complex codelet FFT to input/output array A.

INPUT PARAMETERS:
    A           -   array, must be large enough for plan to work
    Offs        -   offset of the subarray to process
    OperandsCnt -   operands count (see description of FastTransformPlan)
    OperandSize -   operand size (see description of FastTransformPlan)
    MicrovectorSize-microvector size, must be 2
    
OUTPUT PARAMETERS:
    A           -   transformed array

  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftapplycomplexcodeletfft(/* Real    */ ae_vector* a,
     ae_int_t offs,
     ae_int_t operandscnt,
     ae_int_t operandsize,
     ae_int_t microvectorsize,
     ae_state *_state)
{
    ae_int_t opidx;
    ae_int_t n;
    ae_int_t aoffset;
    double a0x;
    double a0y;
    double a1x;
    double a1y;
    double a2x;
    double a2y;
    double a3x;
    double a3y;
    double a4x;
    double a4y;
    double a5x;
    double a5y;
    double v0;
    double v1;
    double v2;
    double v3;
    double t1x;
    double t1y;
    double t2x;
    double t2y;
    double t3x;
    double t3y;
    double t4x;
    double t4y;
    double t5x;
    double t5y;
    double m1x;
    double m1y;
    double m2x;
    double m2y;
    double m3x;
    double m3y;
    double m4x;
    double m4y;
    double m5x;
    double m5y;
    double s1x;
    double s1y;
    double s2x;
    double s2y;
    double s3x;
    double s3y;
    double s4x;
    double s4y;
    double s5x;
    double s5y;
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;
    double v;


    ae_assert(operandscnt>=1, "FTApplyComplexCodeletFFT: OperandsCnt<1", _state);
    ae_assert(operandsize>=1, "FTApplyComplexCodeletFFT: OperandSize<1", _state);
    ae_assert(microvectorsize==2, "FTApplyComplexCodeletFFT: MicrovectorSize<>2", _state);
    n = operandsize;
    
    /*
     * Hard-coded transforms for different N's
     */
    ae_assert(n<=ftbase_maxradix, "FTApplyComplexCodeletFFT: N>MaxRadix", _state);
    if( n==2 )
    {
        for(opidx=0; opidx<=operandscnt-1; opidx++)
        {
            aoffset = offs+opidx*operandsize*2;
            a0x = a->ptr.p_double[aoffset+0];
            a0y = a->ptr.p_double[aoffset+1];
            a1x = a->ptr.p_double[aoffset+2];
            a1y = a->ptr.p_double[aoffset+3];
            v0 = a0x+a1x;
            v1 = a0y+a1y;
            v2 = a0x-a1x;
            v3 = a0y-a1y;
            a->ptr.p_double[aoffset+0] = v0;
            a->ptr.p_double[aoffset+1] = v1;
            a->ptr.p_double[aoffset+2] = v2;
            a->ptr.p_double[aoffset+3] = v3;
        }
        return;
    }
    if( n==3 )
    {
        c1 = ae_cos(2*ae_pi/3, _state)-1;
        c2 = ae_sin(2*ae_pi/3, _state);
        for(opidx=0; opidx<=operandscnt-1; opidx++)
        {
            aoffset = offs+opidx*operandsize*2;
            a0x = a->ptr.p_double[aoffset+0];
            a0y = a->ptr.p_double[aoffset+1];
            a1x = a->ptr.p_double[aoffset+2];
            a1y = a->ptr.p_double[aoffset+3];
            a2x = a->ptr.p_double[aoffset+4];
            a2y = a->ptr.p_double[aoffset+5];
            t1x = a1x+a2x;
            t1y = a1y+a2y;
            a0x = a0x+t1x;
            a0y = a0y+t1y;
            m1x = c1*t1x;
            m1y = c1*t1y;
            m2x = c2*(a1y-a2y);
            m2y = c2*(a2x-a1x);
            s1x = a0x+m1x;
            s1y = a0y+m1y;
            a1x = s1x+m2x;
            a1y = s1y+m2y;
            a2x = s1x-m2x;
            a2y = s1y-m2y;
            a->ptr.p_double[aoffset+0] = a0x;
            a->ptr.p_double[aoffset+1] = a0y;
            a->ptr.p_double[aoffset+2] = a1x;
            a->ptr.p_double[aoffset+3] = a1y;
            a->ptr.p_double[aoffset+4] = a2x;
            a->ptr.p_double[aoffset+5] = a2y;
        }
        return;
    }
    if( n==4 )
    {
        for(opidx=0; opidx<=operandscnt-1; opidx++)
        {
            aoffset = offs+opidx*operandsize*2;
            a0x = a->ptr.p_double[aoffset+0];
            a0y = a->ptr.p_double[aoffset+1];
            a1x = a->ptr.p_double[aoffset+2];
            a1y = a->ptr.p_double[aoffset+3];
            a2x = a->ptr.p_double[aoffset+4];
            a2y = a->ptr.p_double[aoffset+5];
            a3x = a->ptr.p_double[aoffset+6];
            a3y = a->ptr.p_double[aoffset+7];
            t1x = a0x+a2x;
            t1y = a0y+a2y;
            t2x = a1x+a3x;
            t2y = a1y+a3y;
            m2x = a0x-a2x;
            m2y = a0y-a2y;
            m3x = a1y-a3y;
            m3y = a3x-a1x;
            a->ptr.p_double[aoffset+0] = t1x+t2x;
            a->ptr.p_double[aoffset+1] = t1y+t2y;
            a->ptr.p_double[aoffset+4] = t1x-t2x;
            a->ptr.p_double[aoffset+5] = t1y-t2y;
            a->ptr.p_double[aoffset+2] = m2x+m3x;
            a->ptr.p_double[aoffset+3] = m2y+m3y;
            a->ptr.p_double[aoffset+6] = m2x-m3x;
            a->ptr.p_double[aoffset+7] = m2y-m3y;
        }
        return;
    }
    if( n==5 )
    {
        v = 2*ae_pi/5;
        c1 = (ae_cos(v, _state)+ae_cos(2*v, _state))/2-1;
        c2 = (ae_cos(v, _state)-ae_cos(2*v, _state))/2;
        c3 = -ae_sin(v, _state);
        c4 = -(ae_sin(v, _state)+ae_sin(2*v, _state));
        c5 = ae_sin(v, _state)-ae_sin(2*v, _state);
        for(opidx=0; opidx<=operandscnt-1; opidx++)
        {
            aoffset = offs+opidx*operandsize*2;
            t1x = a->ptr.p_double[aoffset+2]+a->ptr.p_double[aoffset+8];
            t1y = a->ptr.p_double[aoffset+3]+a->ptr.p_double[aoffset+9];
            t2x = a->ptr.p_double[aoffset+4]+a->ptr.p_double[aoffset+6];
            t2y = a->ptr.p_double[aoffset+5]+a->ptr.p_double[aoffset+7];
            t3x = a->ptr.p_double[aoffset+2]-a->ptr.p_double[aoffset+8];
            t3y = a->ptr.p_double[aoffset+3]-a->ptr.p_double[aoffset+9];
            t4x = a->ptr.p_double[aoffset+6]-a->ptr.p_double[aoffset+4];
            t4y = a->ptr.p_double[aoffset+7]-a->ptr.p_double[aoffset+5];
            t5x = t1x+t2x;
            t5y = t1y+t2y;
            a->ptr.p_double[aoffset+0] = a->ptr.p_double[aoffset+0]+t5x;
            a->ptr.p_double[aoffset+1] = a->ptr.p_double[aoffset+1]+t5y;
            m1x = c1*t5x;
            m1y = c1*t5y;
            m2x = c2*(t1x-t2x);
            m2y = c2*(t1y-t2y);
            m3x = -c3*(t3y+t4y);
            m3y = c3*(t3x+t4x);
            m4x = -c4*t4y;
            m4y = c4*t4x;
            m5x = -c5*t3y;
            m5y = c5*t3x;
            s3x = m3x-m4x;
            s3y = m3y-m4y;
            s5x = m3x+m5x;
            s5y = m3y+m5y;
            s1x = a->ptr.p_double[aoffset+0]+m1x;
            s1y = a->ptr.p_double[aoffset+1]+m1y;
            s2x = s1x+m2x;
            s2y = s1y+m2y;
            s4x = s1x-m2x;
            s4y = s1y-m2y;
            a->ptr.p_double[aoffset+2] = s2x+s3x;
            a->ptr.p_double[aoffset+3] = s2y+s3y;
            a->ptr.p_double[aoffset+4] = s4x+s5x;
            a->ptr.p_double[aoffset+5] = s4y+s5y;
            a->ptr.p_double[aoffset+6] = s4x-s5x;
            a->ptr.p_double[aoffset+7] = s4y-s5y;
            a->ptr.p_double[aoffset+8] = s2x-s3x;
            a->ptr.p_double[aoffset+9] = s2y-s3y;
        }
        return;
    }
    if( n==6 )
    {
        c1 = ae_cos(2*ae_pi/3, _state)-1;
        c2 = ae_sin(2*ae_pi/3, _state);
        c3 = ae_cos(-ae_pi/3, _state);
        c4 = ae_sin(-ae_pi/3, _state);
        for(opidx=0; opidx<=operandscnt-1; opidx++)
        {
            aoffset = offs+opidx*operandsize*2;
            a0x = a->ptr.p_double[aoffset+0];
            a0y = a->ptr.p_double[aoffset+1];
            a1x = a->ptr.p_double[aoffset+2];
            a1y = a->ptr.p_double[aoffset+3];
            a2x = a->ptr.p_double[aoffset+4];
            a2y = a->ptr.p_double[aoffset+5];
            a3x = a->ptr.p_double[aoffset+6];
            a3y = a->ptr.p_double[aoffset+7];
            a4x = a->ptr.p_double[aoffset+8];
            a4y = a->ptr.p_double[aoffset+9];
            a5x = a->ptr.p_double[aoffset+10];
            a5y = a->ptr.p_double[aoffset+11];
            v0 = a0x;
            v1 = a0y;
            a0x = a0x+a3x;
            a0y = a0y+a3y;
            a3x = v0-a3x;
            a3y = v1-a3y;
            v0 = a1x;
            v1 = a1y;
            a1x = a1x+a4x;
            a1y = a1y+a4y;
            a4x = v0-a4x;
            a4y = v1-a4y;
            v0 = a2x;
            v1 = a2y;
            a2x = a2x+a5x;
            a2y = a2y+a5y;
            a5x = v0-a5x;
            a5y = v1-a5y;
            t4x = a4x*c3-a4y*c4;
            t4y = a4x*c4+a4y*c3;
            a4x = t4x;
            a4y = t4y;
            t5x = -a5x*c3-a5y*c4;
            t5y = a5x*c4-a5y*c3;
            a5x = t5x;
            a5y = t5y;
            t1x = a1x+a2x;
            t1y = a1y+a2y;
            a0x = a0x+t1x;
            a0y = a0y+t1y;
            m1x = c1*t1x;
            m1y = c1*t1y;
            m2x = c2*(a1y-a2y);
            m2y = c2*(a2x-a1x);
            s1x = a0x+m1x;
            s1y = a0y+m1y;
            a1x = s1x+m2x;
            a1y = s1y+m2y;
            a2x = s1x-m2x;
            a2y = s1y-m2y;
            t1x = a4x+a5x;
            t1y = a4y+a5y;
            a3x = a3x+t1x;
            a3y = a3y+t1y;
            m1x = c1*t1x;
            m1y = c1*t1y;
            m2x = c2*(a4y-a5y);
            m2y = c2*(a5x-a4x);
            s1x = a3x+m1x;
            s1y = a3y+m1y;
            a4x = s1x+m2x;
            a4y = s1y+m2y;
            a5x = s1x-m2x;
            a5y = s1y-m2y;
            a->ptr.p_double[aoffset+0] = a0x;
            a->ptr.p_double[aoffset+1] = a0y;
            a->ptr.p_double[aoffset+2] = a3x;
            a->ptr.p_double[aoffset+3] = a3y;
            a->ptr.p_double[aoffset+4] = a1x;
            a->ptr.p_double[aoffset+5] = a1y;
            a->ptr.p_double[aoffset+6] = a4x;
            a->ptr.p_double[aoffset+7] = a4y;
            a->ptr.p_double[aoffset+8] = a2x;
            a->ptr.p_double[aoffset+9] = a2y;
            a->ptr.p_double[aoffset+10] = a5x;
            a->ptr.p_double[aoffset+11] = a5y;
        }
        return;
    }
}


/*************************************************************************
This subroutine applies complex "integrated" codelet FFT  to  input/output
array A. "Integrated" codelet differs from "normal" one in following ways:
* it can work with MicrovectorSize>1
* hence, it can be used in Cooley-Tukey FFT without transpositions
* it performs inlined multiplication by twiddle factors of Cooley-Tukey
  FFT with N2=MicrovectorSize/2.

INPUT PARAMETERS:
    A           -   array, must be large enough for plan to work
    Offs        -   offset of the subarray to process
    OperandsCnt -   operands count (see description of FastTransformPlan)
    OperandSize -   operand size (see description of FastTransformPlan)
    MicrovectorSize-microvector size, must be 1
    
OUTPUT PARAMETERS:
    A           -   transformed array

  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftapplycomplexcodelettwfft(/* Real    */ ae_vector* a,
     ae_int_t offs,
     ae_int_t operandscnt,
     ae_int_t operandsize,
     ae_int_t microvectorsize,
     ae_state *_state)
{
    ae_int_t opidx;
    ae_int_t mvidx;
    ae_int_t n;
    ae_int_t m;
    ae_int_t aoffset0;
    ae_int_t aoffset2;
    ae_int_t aoffset4;
    ae_int_t aoffset6;
    ae_int_t aoffset8;
    ae_int_t aoffset10;
    double a0x;
    double a0y;
    double a1x;
    double a1y;
    double a2x;
    double a2y;
    double a3x;
    double a3y;
    double a4x;
    double a4y;
    double a5x;
    double a5y;
    double v0;
    double v1;
    double v2;
    double v3;
    double q0x;
    double q0y;
    double t1x;
    double t1y;
    double t2x;
    double t2y;
    double t3x;
    double t3y;
    double t4x;
    double t4y;
    double t5x;
    double t5y;
    double m1x;
    double m1y;
    double m2x;
    double m2y;
    double m3x;
    double m3y;
    double m4x;
    double m4y;
    double m5x;
    double m5y;
    double s1x;
    double s1y;
    double s2x;
    double s2y;
    double s3x;
    double s3y;
    double s4x;
    double s4y;
    double s5x;
    double s5y;
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;
    double v;
    double tw0;
    double tw1;
    double twx;
    double twxm1;
    double twy;
    double tw2x;
    double tw2y;
    double tw3x;
    double tw3y;
    double tw4x;
    double tw4y;
    double tw5x;
    double tw5y;


    ae_assert(operandscnt>=1, "FTApplyComplexCodeletFFT: OperandsCnt<1", _state);
    ae_assert(operandsize>=1, "FTApplyComplexCodeletFFT: OperandSize<1", _state);
    ae_assert(microvectorsize>=1, "FTApplyComplexCodeletFFT: MicrovectorSize<>1", _state);
    ae_assert(microvectorsize%2==0, "FTApplyComplexCodeletFFT: MicrovectorSize is not even", _state);
    n = operandsize;
    m = microvectorsize/2;
    
    /*
     * Hard-coded transforms for different N's
     */
    ae_assert(n<=ftbase_maxradix, "FTApplyComplexCodeletTwFFT: N>MaxRadix", _state);
    if( n==2 )
    {
        v = -2*ae_pi/(n*m);
        tw0 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
        tw1 = ae_sin(v, _state);
        for(opidx=0; opidx<=operandscnt-1; opidx++)
        {
            aoffset0 = offs+opidx*operandsize*microvectorsize;
            aoffset2 = aoffset0+microvectorsize;
            twxm1 = 0.0;
            twy = 0.0;
            for(mvidx=0; mvidx<=m-1; mvidx++)
            {
                a0x = a->ptr.p_double[aoffset0];
                a0y = a->ptr.p_double[aoffset0+1];
                a1x = a->ptr.p_double[aoffset2];
                a1y = a->ptr.p_double[aoffset2+1];
                v0 = a0x+a1x;
                v1 = a0y+a1y;
                v2 = a0x-a1x;
                v3 = a0y-a1y;
                a->ptr.p_double[aoffset0] = v0;
                a->ptr.p_double[aoffset0+1] = v1;
                a->ptr.p_double[aoffset2] = v2*(1+twxm1)-v3*twy;
                a->ptr.p_double[aoffset2+1] = v3*(1+twxm1)+v2*twy;
                aoffset0 = aoffset0+2;
                aoffset2 = aoffset2+2;
                if( (mvidx+1)%ftbase_updatetw==0 )
                {
                    v = -2*ae_pi*(mvidx+1)/(n*m);
                    twxm1 = ae_sin(0.5*v, _state);
                    twxm1 = -2*twxm1*twxm1;
                    twy = ae_sin(v, _state);
                }
                else
                {
                    v = twxm1+tw0+twxm1*tw0-twy*tw1;
                    twy = twy+tw1+twxm1*tw1+twy*tw0;
                    twxm1 = v;
                }
            }
        }
        return;
    }
    if( n==3 )
    {
        v = -2*ae_pi/(n*m);
        tw0 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
        tw1 = ae_sin(v, _state);
        c1 = ae_cos(2*ae_pi/3, _state)-1;
        c2 = ae_sin(2*ae_pi/3, _state);
        for(opidx=0; opidx<=operandscnt-1; opidx++)
        {
            aoffset0 = offs+opidx*operandsize*microvectorsize;
            aoffset2 = aoffset0+microvectorsize;
            aoffset4 = aoffset2+microvectorsize;
            twx = 1.0;
            twxm1 = 0.0;
            twy = 0.0;
            for(mvidx=0; mvidx<=m-1; mvidx++)
            {
                a0x = a->ptr.p_double[aoffset0];
                a0y = a->ptr.p_double[aoffset0+1];
                a1x = a->ptr.p_double[aoffset2];
                a1y = a->ptr.p_double[aoffset2+1];
                a2x = a->ptr.p_double[aoffset4];
                a2y = a->ptr.p_double[aoffset4+1];
                t1x = a1x+a2x;
                t1y = a1y+a2y;
                a0x = a0x+t1x;
                a0y = a0y+t1y;
                m1x = c1*t1x;
                m1y = c1*t1y;
                m2x = c2*(a1y-a2y);
                m2y = c2*(a2x-a1x);
                s1x = a0x+m1x;
                s1y = a0y+m1y;
                a1x = s1x+m2x;
                a1y = s1y+m2y;
                a2x = s1x-m2x;
                a2y = s1y-m2y;
                tw2x = twx*twx-twy*twy;
                tw2y = 2*twx*twy;
                a->ptr.p_double[aoffset0] = a0x;
                a->ptr.p_double[aoffset0+1] = a0y;
                a->ptr.p_double[aoffset2] = a1x*twx-a1y*twy;
                a->ptr.p_double[aoffset2+1] = a1y*twx+a1x*twy;
                a->ptr.p_double[aoffset4] = a2x*tw2x-a2y*tw2y;
                a->ptr.p_double[aoffset4+1] = a2y*tw2x+a2x*tw2y;
                aoffset0 = aoffset0+2;
                aoffset2 = aoffset2+2;
                aoffset4 = aoffset4+2;
                if( (mvidx+1)%ftbase_updatetw==0 )
                {
                    v = -2*ae_pi*(mvidx+1)/(n*m);
                    twxm1 = ae_sin(0.5*v, _state);
                    twxm1 = -2*twxm1*twxm1;
                    twy = ae_sin(v, _state);
                    twx = twxm1+1;
                }
                else
                {
                    v = twxm1+tw0+twxm1*tw0-twy*tw1;
                    twy = twy+tw1+twxm1*tw1+twy*tw0;
                    twxm1 = v;
                    twx = v+1;
                }
            }
        }
        return;
    }
    if( n==4 )
    {
        v = -2*ae_pi/(n*m);
        tw0 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
        tw1 = ae_sin(v, _state);
        for(opidx=0; opidx<=operandscnt-1; opidx++)
        {
            aoffset0 = offs+opidx*operandsize*microvectorsize;
            aoffset2 = aoffset0+microvectorsize;
            aoffset4 = aoffset2+microvectorsize;
            aoffset6 = aoffset4+microvectorsize;
            twx = 1.0;
            twxm1 = 0.0;
            twy = 0.0;
            for(mvidx=0; mvidx<=m-1; mvidx++)
            {
                a0x = a->ptr.p_double[aoffset0];
                a0y = a->ptr.p_double[aoffset0+1];
                a1x = a->ptr.p_double[aoffset2];
                a1y = a->ptr.p_double[aoffset2+1];
                a2x = a->ptr.p_double[aoffset4];
                a2y = a->ptr.p_double[aoffset4+1];
                a3x = a->ptr.p_double[aoffset6];
                a3y = a->ptr.p_double[aoffset6+1];
                t1x = a0x+a2x;
                t1y = a0y+a2y;
                t2x = a1x+a3x;
                t2y = a1y+a3y;
                m2x = a0x-a2x;
                m2y = a0y-a2y;
                m3x = a1y-a3y;
                m3y = a3x-a1x;
                tw2x = twx*twx-twy*twy;
                tw2y = 2*twx*twy;
                tw3x = twx*tw2x-twy*tw2y;
                tw3y = twx*tw2y+twy*tw2x;
                a1x = m2x+m3x;
                a1y = m2y+m3y;
                a2x = t1x-t2x;
                a2y = t1y-t2y;
                a3x = m2x-m3x;
                a3y = m2y-m3y;
                a->ptr.p_double[aoffset0] = t1x+t2x;
                a->ptr.p_double[aoffset0+1] = t1y+t2y;
                a->ptr.p_double[aoffset2] = a1x*twx-a1y*twy;
                a->ptr.p_double[aoffset2+1] = a1y*twx+a1x*twy;
                a->ptr.p_double[aoffset4] = a2x*tw2x-a2y*tw2y;
                a->ptr.p_double[aoffset4+1] = a2y*tw2x+a2x*tw2y;
                a->ptr.p_double[aoffset6] = a3x*tw3x-a3y*tw3y;
                a->ptr.p_double[aoffset6+1] = a3y*tw3x+a3x*tw3y;
                aoffset0 = aoffset0+2;
                aoffset2 = aoffset2+2;
                aoffset4 = aoffset4+2;
                aoffset6 = aoffset6+2;
                if( (mvidx+1)%ftbase_updatetw==0 )
                {
                    v = -2*ae_pi*(mvidx+1)/(n*m);
                    twxm1 = ae_sin(0.5*v, _state);
                    twxm1 = -2*twxm1*twxm1;
                    twy = ae_sin(v, _state);
                    twx = twxm1+1;
                }
                else
                {
                    v = twxm1+tw0+twxm1*tw0-twy*tw1;
                    twy = twy+tw1+twxm1*tw1+twy*tw0;
                    twxm1 = v;
                    twx = v+1;
                }
            }
        }
        return;
    }
    if( n==5 )
    {
        v = -2*ae_pi/(n*m);
        tw0 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
        tw1 = ae_sin(v, _state);
        v = 2*ae_pi/5;
        c1 = (ae_cos(v, _state)+ae_cos(2*v, _state))/2-1;
        c2 = (ae_cos(v, _state)-ae_cos(2*v, _state))/2;
        c3 = -ae_sin(v, _state);
        c4 = -(ae_sin(v, _state)+ae_sin(2*v, _state));
        c5 = ae_sin(v, _state)-ae_sin(2*v, _state);
        for(opidx=0; opidx<=operandscnt-1; opidx++)
        {
            aoffset0 = offs+opidx*operandsize*microvectorsize;
            aoffset2 = aoffset0+microvectorsize;
            aoffset4 = aoffset2+microvectorsize;
            aoffset6 = aoffset4+microvectorsize;
            aoffset8 = aoffset6+microvectorsize;
            twx = 1.0;
            twxm1 = 0.0;
            twy = 0.0;
            for(mvidx=0; mvidx<=m-1; mvidx++)
            {
                a0x = a->ptr.p_double[aoffset0];
                a0y = a->ptr.p_double[aoffset0+1];
                a1x = a->ptr.p_double[aoffset2];
                a1y = a->ptr.p_double[aoffset2+1];
                a2x = a->ptr.p_double[aoffset4];
                a2y = a->ptr.p_double[aoffset4+1];
                a3x = a->ptr.p_double[aoffset6];
                a3y = a->ptr.p_double[aoffset6+1];
                a4x = a->ptr.p_double[aoffset8];
                a4y = a->ptr.p_double[aoffset8+1];
                t1x = a1x+a4x;
                t1y = a1y+a4y;
                t2x = a2x+a3x;
                t2y = a2y+a3y;
                t3x = a1x-a4x;
                t3y = a1y-a4y;
                t4x = a3x-a2x;
                t4y = a3y-a2y;
                t5x = t1x+t2x;
                t5y = t1y+t2y;
                q0x = a0x+t5x;
                q0y = a0y+t5y;
                m1x = c1*t5x;
                m1y = c1*t5y;
                m2x = c2*(t1x-t2x);
                m2y = c2*(t1y-t2y);
                m3x = -c3*(t3y+t4y);
                m3y = c3*(t3x+t4x);
                m4x = -c4*t4y;
                m4y = c4*t4x;
                m5x = -c5*t3y;
                m5y = c5*t3x;
                s3x = m3x-m4x;
                s3y = m3y-m4y;
                s5x = m3x+m5x;
                s5y = m3y+m5y;
                s1x = q0x+m1x;
                s1y = q0y+m1y;
                s2x = s1x+m2x;
                s2y = s1y+m2y;
                s4x = s1x-m2x;
                s4y = s1y-m2y;
                tw2x = twx*twx-twy*twy;
                tw2y = 2*twx*twy;
                tw3x = twx*tw2x-twy*tw2y;
                tw3y = twx*tw2y+twy*tw2x;
                tw4x = tw2x*tw2x-tw2y*tw2y;
                tw4y = tw2x*tw2y+tw2y*tw2x;
                a1x = s2x+s3x;
                a1y = s2y+s3y;
                a2x = s4x+s5x;
                a2y = s4y+s5y;
                a3x = s4x-s5x;
                a3y = s4y-s5y;
                a4x = s2x-s3x;
                a4y = s2y-s3y;
                a->ptr.p_double[aoffset0] = q0x;
                a->ptr.p_double[aoffset0+1] = q0y;
                a->ptr.p_double[aoffset2] = a1x*twx-a1y*twy;
                a->ptr.p_double[aoffset2+1] = a1x*twy+a1y*twx;
                a->ptr.p_double[aoffset4] = a2x*tw2x-a2y*tw2y;
                a->ptr.p_double[aoffset4+1] = a2x*tw2y+a2y*tw2x;
                a->ptr.p_double[aoffset6] = a3x*tw3x-a3y*tw3y;
                a->ptr.p_double[aoffset6+1] = a3x*tw3y+a3y*tw3x;
                a->ptr.p_double[aoffset8] = a4x*tw4x-a4y*tw4y;
                a->ptr.p_double[aoffset8+1] = a4x*tw4y+a4y*tw4x;
                aoffset0 = aoffset0+2;
                aoffset2 = aoffset2+2;
                aoffset4 = aoffset4+2;
                aoffset6 = aoffset6+2;
                aoffset8 = aoffset8+2;
                if( (mvidx+1)%ftbase_updatetw==0 )
                {
                    v = -2*ae_pi*(mvidx+1)/(n*m);
                    twxm1 = ae_sin(0.5*v, _state);
                    twxm1 = -2*twxm1*twxm1;
                    twy = ae_sin(v, _state);
                    twx = twxm1+1;
                }
                else
                {
                    v = twxm1+tw0+twxm1*tw0-twy*tw1;
                    twy = twy+tw1+twxm1*tw1+twy*tw0;
                    twxm1 = v;
                    twx = v+1;
                }
            }
        }
        return;
    }
    if( n==6 )
    {
        c1 = ae_cos(2*ae_pi/3, _state)-1;
        c2 = ae_sin(2*ae_pi/3, _state);
        c3 = ae_cos(-ae_pi/3, _state);
        c4 = ae_sin(-ae_pi/3, _state);
        v = -2*ae_pi/(n*m);
        tw0 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
        tw1 = ae_sin(v, _state);
        for(opidx=0; opidx<=operandscnt-1; opidx++)
        {
            aoffset0 = offs+opidx*operandsize*microvectorsize;
            aoffset2 = aoffset0+microvectorsize;
            aoffset4 = aoffset2+microvectorsize;
            aoffset6 = aoffset4+microvectorsize;
            aoffset8 = aoffset6+microvectorsize;
            aoffset10 = aoffset8+microvectorsize;
            twx = 1.0;
            twxm1 = 0.0;
            twy = 0.0;
            for(mvidx=0; mvidx<=m-1; mvidx++)
            {
                a0x = a->ptr.p_double[aoffset0+0];
                a0y = a->ptr.p_double[aoffset0+1];
                a1x = a->ptr.p_double[aoffset2+0];
                a1y = a->ptr.p_double[aoffset2+1];
                a2x = a->ptr.p_double[aoffset4+0];
                a2y = a->ptr.p_double[aoffset4+1];
                a3x = a->ptr.p_double[aoffset6+0];
                a3y = a->ptr.p_double[aoffset6+1];
                a4x = a->ptr.p_double[aoffset8+0];
                a4y = a->ptr.p_double[aoffset8+1];
                a5x = a->ptr.p_double[aoffset10+0];
                a5y = a->ptr.p_double[aoffset10+1];
                v0 = a0x;
                v1 = a0y;
                a0x = a0x+a3x;
                a0y = a0y+a3y;
                a3x = v0-a3x;
                a3y = v1-a3y;
                v0 = a1x;
                v1 = a1y;
                a1x = a1x+a4x;
                a1y = a1y+a4y;
                a4x = v0-a4x;
                a4y = v1-a4y;
                v0 = a2x;
                v1 = a2y;
                a2x = a2x+a5x;
                a2y = a2y+a5y;
                a5x = v0-a5x;
                a5y = v1-a5y;
                t4x = a4x*c3-a4y*c4;
                t4y = a4x*c4+a4y*c3;
                a4x = t4x;
                a4y = t4y;
                t5x = -a5x*c3-a5y*c4;
                t5y = a5x*c4-a5y*c3;
                a5x = t5x;
                a5y = t5y;
                t1x = a1x+a2x;
                t1y = a1y+a2y;
                a0x = a0x+t1x;
                a0y = a0y+t1y;
                m1x = c1*t1x;
                m1y = c1*t1y;
                m2x = c2*(a1y-a2y);
                m2y = c2*(a2x-a1x);
                s1x = a0x+m1x;
                s1y = a0y+m1y;
                a1x = s1x+m2x;
                a1y = s1y+m2y;
                a2x = s1x-m2x;
                a2y = s1y-m2y;
                t1x = a4x+a5x;
                t1y = a4y+a5y;
                a3x = a3x+t1x;
                a3y = a3y+t1y;
                m1x = c1*t1x;
                m1y = c1*t1y;
                m2x = c2*(a4y-a5y);
                m2y = c2*(a5x-a4x);
                s1x = a3x+m1x;
                s1y = a3y+m1y;
                a4x = s1x+m2x;
                a4y = s1y+m2y;
                a5x = s1x-m2x;
                a5y = s1y-m2y;
                tw2x = twx*twx-twy*twy;
                tw2y = 2*twx*twy;
                tw3x = twx*tw2x-twy*tw2y;
                tw3y = twx*tw2y+twy*tw2x;
                tw4x = tw2x*tw2x-tw2y*tw2y;
                tw4y = 2*tw2x*tw2y;
                tw5x = tw3x*tw2x-tw3y*tw2y;
                tw5y = tw3x*tw2y+tw3y*tw2x;
                a->ptr.p_double[aoffset0+0] = a0x;
                a->ptr.p_double[aoffset0+1] = a0y;
                a->ptr.p_double[aoffset2+0] = a3x*twx-a3y*twy;
                a->ptr.p_double[aoffset2+1] = a3y*twx+a3x*twy;
                a->ptr.p_double[aoffset4+0] = a1x*tw2x-a1y*tw2y;
                a->ptr.p_double[aoffset4+1] = a1y*tw2x+a1x*tw2y;
                a->ptr.p_double[aoffset6+0] = a4x*tw3x-a4y*tw3y;
                a->ptr.p_double[aoffset6+1] = a4y*tw3x+a4x*tw3y;
                a->ptr.p_double[aoffset8+0] = a2x*tw4x-a2y*tw4y;
                a->ptr.p_double[aoffset8+1] = a2y*tw4x+a2x*tw4y;
                a->ptr.p_double[aoffset10+0] = a5x*tw5x-a5y*tw5y;
                a->ptr.p_double[aoffset10+1] = a5y*tw5x+a5x*tw5y;
                aoffset0 = aoffset0+2;
                aoffset2 = aoffset2+2;
                aoffset4 = aoffset4+2;
                aoffset6 = aoffset6+2;
                aoffset8 = aoffset8+2;
                aoffset10 = aoffset10+2;
                if( (mvidx+1)%ftbase_updatetw==0 )
                {
                    v = -2*ae_pi*(mvidx+1)/(n*m);
                    twxm1 = ae_sin(0.5*v, _state);
                    twxm1 = -2*twxm1*twxm1;
                    twy = ae_sin(v, _state);
                    twx = twxm1+1;
                }
                else
                {
                    v = twxm1+tw0+twxm1*tw0-twy*tw1;
                    twy = twy+tw1+twxm1*tw1+twy*tw0;
                    twxm1 = v;
                    twx = v+1;
                }
            }
        }
        return;
    }
}


/*************************************************************************
This subroutine precomputes data for complex Bluestein's  FFT  and  writes
them to array PrecR[] at specified offset. It  is  responsibility  of  the
caller to make sure that PrecR[] is large enough.

INPUT PARAMETERS:
    N           -   original size of the transform
    M           -   size of the "padded" Bluestein's transform
    PrecR       -   preallocated array
    Offs        -   offset
    
OUTPUT PARAMETERS:
    PrecR       -   data at Offs:Offs+4*M-1 are modified:
                    * PrecR[Offs:Offs+2*M-1] stores Z[k]=exp(i*pi*k^2/N)
                    * PrecR[Offs+2*M:Offs+4*M-1] stores FFT of the Z
                    Other parts of PrecR are unchanged.
                    
NOTE: this function performs internal M-point FFT. It allocates temporary
      plan which is destroyed after leaving this function.

  -- ALGLIB --
     Copyright 08.05.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftprecomputebluesteinsfft(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* precr,
     ae_int_t offs,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    double bx;
    double by;
    fasttransformplan plan;

    ae_frame_make(_state, &_frame_block);
    _fasttransformplan_init(&plan, _state);

    
    /*
     * Fill first half of PrecR with b[k] = exp(i*pi*k^2/N)
     */
    for(i=0; i<=2*m-1; i++)
    {
        precr->ptr.p_double[offs+i] = (double)(0);
    }
    for(i=0; i<=n-1; i++)
    {
        bx = ae_cos(ae_pi/n*i*i, _state);
        by = ae_sin(ae_pi/n*i*i, _state);
        precr->ptr.p_double[offs+2*i+0] = bx;
        precr->ptr.p_double[offs+2*i+1] = by;
        precr->ptr.p_double[offs+2*((m-i)%m)+0] = bx;
        precr->ptr.p_double[offs+2*((m-i)%m)+1] = by;
    }
    
    /*
     * Precomputed FFT
     */
    ftcomplexfftplan(m, 1, &plan, _state);
    for(i=0; i<=2*m-1; i++)
    {
        precr->ptr.p_double[offs+2*m+i] = precr->ptr.p_double[offs+i];
    }
    ftbase_ftapplysubplan(&plan, 0, precr, offs+2*m, 0, &plan.buffer, 1, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
This subroutine applies complex Bluestein's FFT to input/output array A.

INPUT PARAMETERS:
    Plan        -   transformation plan
    A           -   array, must be large enough for plan to work
    ABase       -   base offset in array A, this value points to start of
                    subarray whose length is equal to length of the plan
    AOffset     -   offset with respect to ABase, 0<=AOffset<PlanLength.
                    This is an offset within large PlanLength-subarray of
                    the chunk to process.
    OperandsCnt -   number of repeated operands (length N each)
    N           -   original data length (measured in complex numbers)
    M           -   padded data length (measured in complex numbers)
    PrecOffs    -   offset of the precomputed data for the plan
    SubPlan     -   position of the length-M FFT subplan which is used by
                    transformation
    BufA        -   temporary buffer, at least 2*M elements
    BufB        -   temporary buffer, at least 2*M elements
    BufC        -   temporary buffer, at least 2*M elements
    BufD        -   temporary buffer, at least 2*M elements
    
OUTPUT PARAMETERS:
    A           -   transformed array

  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftbluesteinsfft(fasttransformplan* plan,
     /* Real    */ ae_vector* a,
     ae_int_t abase,
     ae_int_t aoffset,
     ae_int_t operandscnt,
     ae_int_t n,
     ae_int_t m,
     ae_int_t precoffs,
     ae_int_t subplan,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     /* Real    */ ae_vector* bufc,
     /* Real    */ ae_vector* bufd,
     ae_state *_state)
{
    ae_int_t op;
    ae_int_t i;
    double x;
    double y;
    double bx;
    double by;
    double ax;
    double ay;
    double rx;
    double ry;
    ae_int_t p0;
    ae_int_t p1;
    ae_int_t p2;


    for(op=0; op<=operandscnt-1; op++)
    {
        
        /*
         * Multiply A by conj(Z), store to buffer.
         * Pad A by zeros.
         *
         * NOTE: Z[k]=exp(i*pi*k^2/N)
         */
        p0 = abase+aoffset+op*2*n;
        p1 = precoffs;
        for(i=0; i<=n-1; i++)
        {
            x = a->ptr.p_double[p0+0];
            y = a->ptr.p_double[p0+1];
            bx = plan->precr.ptr.p_double[p1+0];
            by = -plan->precr.ptr.p_double[p1+1];
            bufa->ptr.p_double[2*i+0] = x*bx-y*by;
            bufa->ptr.p_double[2*i+1] = x*by+y*bx;
            p0 = p0+2;
            p1 = p1+2;
        }
        for(i=2*n; i<=2*m-1; i++)
        {
            bufa->ptr.p_double[i] = (double)(0);
        }
        
        /*
         * Perform convolution of A and Z (using precomputed
         * FFT of Z stored in Plan structure).
         */
        ftbase_ftapplysubplan(plan, subplan, bufa, 0, 0, bufc, 1, _state);
        p0 = 0;
        p1 = precoffs+2*m;
        for(i=0; i<=m-1; i++)
        {
            ax = bufa->ptr.p_double[p0+0];
            ay = bufa->ptr.p_double[p0+1];
            bx = plan->precr.ptr.p_double[p1+0];
            by = plan->precr.ptr.p_double[p1+1];
            bufa->ptr.p_double[p0+0] = ax*bx-ay*by;
            bufa->ptr.p_double[p0+1] = -(ax*by+ay*bx);
            p0 = p0+2;
            p1 = p1+2;
        }
        ftbase_ftapplysubplan(plan, subplan, bufa, 0, 0, bufc, 1, _state);
        
        /*
         * Post processing:
         *     A:=conj(Z)*conj(A)/M
         * Here conj(A)/M corresponds to last stage of inverse DFT,
         * and conj(Z) comes from Bluestein's FFT algorithm.
         */
        p0 = precoffs;
        p1 = 0;
        p2 = abase+aoffset+op*2*n;
        for(i=0; i<=n-1; i++)
        {
            bx = plan->precr.ptr.p_double[p0+0];
            by = plan->precr.ptr.p_double[p0+1];
            rx = bufa->ptr.p_double[p1+0]/m;
            ry = -bufa->ptr.p_double[p1+1]/m;
            a->ptr.p_double[p2+0] = rx*bx-ry*(-by);
            a->ptr.p_double[p2+1] = rx*(-by)+ry*bx;
            p0 = p0+2;
            p1 = p1+2;
            p2 = p2+2;
        }
    }
}


/*************************************************************************
This subroutine precomputes data for complex Rader's FFT and  writes  them
to array PrecR[] at specified offset. It  is  responsibility of the caller
to make sure that PrecR[] is large enough.

INPUT PARAMETERS:
    N           -   original size of the transform (before reduction to N-1)
    RQ          -   primitive root modulo N
    RIQ         -   inverse of primitive root modulo N
    PrecR       -   preallocated array
    Offs        -   offset
    
OUTPUT PARAMETERS:
    PrecR       -   data at Offs:Offs+2*(N-1)-1 store FFT of Rader's factors,
                    other parts of PrecR are unchanged.
                    
NOTE: this function performs internal (N-1)-point FFT. It allocates temporary
      plan which is destroyed after leaving this function.

  -- ALGLIB --
     Copyright 08.05.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftprecomputeradersfft(ae_int_t n,
     ae_int_t rq,
     ae_int_t riq,
     /* Real    */ ae_vector* precr,
     ae_int_t offs,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t q;
    fasttransformplan plan;
    ae_int_t kiq;
    double v;

    ae_frame_make(_state, &_frame_block);
    _fasttransformplan_init(&plan, _state);

    
    /*
     * Fill PrecR with Rader factors, perform FFT
     */
    kiq = 1;
    for(q=0; q<=n-2; q++)
    {
        v = -2*ae_pi*kiq/n;
        precr->ptr.p_double[offs+2*q+0] = ae_cos(v, _state);
        precr->ptr.p_double[offs+2*q+1] = ae_sin(v, _state);
        kiq = kiq*riq%n;
    }
    ftcomplexfftplan(n-1, 1, &plan, _state);
    ftbase_ftapplysubplan(&plan, 0, precr, offs, 0, &plan.buffer, 1, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
This subroutine applies complex Rader's FFT to input/output array A.

INPUT PARAMETERS:
    A           -   array, must be large enough for plan to work
    ABase       -   base offset in array A, this value points to start of
                    subarray whose length is equal to length of the plan
    AOffset     -   offset with respect to ABase, 0<=AOffset<PlanLength.
                    This is an offset within large PlanLength-subarray of
                    the chunk to process.
    OperandsCnt -   number of repeated operands (length N each)
    N           -   original data length (measured in complex numbers)
    SubPlan     -   position of the (N-1)-point FFT subplan which is used
                    by transformation
    RQ          -   primitive root modulo N
    RIQ         -   inverse of primitive root modulo N
    PrecOffs    -   offset of the precomputed data for the plan
    Buf         -   temporary array
    
OUTPUT PARAMETERS:
    A           -   transformed array

  -- ALGLIB --
     Copyright 05.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftradersfft(fasttransformplan* plan,
     /* Real    */ ae_vector* a,
     ae_int_t abase,
     ae_int_t aoffset,
     ae_int_t operandscnt,
     ae_int_t n,
     ae_int_t subplan,
     ae_int_t rq,
     ae_int_t riq,
     ae_int_t precoffs,
     /* Real    */ ae_vector* buf,
     ae_state *_state)
{
    ae_int_t opidx;
    ae_int_t i;
    ae_int_t q;
    ae_int_t kq;
    ae_int_t kiq;
    double x0;
    double y0;
    ae_int_t p0;
    ae_int_t p1;
    double ax;
    double ay;
    double bx;
    double by;
    double rx;
    double ry;


    ae_assert(operandscnt>=1, "FTApplyComplexRefFFT: OperandsCnt<1", _state);
    
    /*
     * Process operands
     */
    for(opidx=0; opidx<=operandscnt-1; opidx++)
    {
        
        /*
         * fill QA
         */
        kq = 1;
        p0 = abase+aoffset+opidx*n*2;
        p1 = aoffset+opidx*n*2;
        rx = a->ptr.p_double[p0+0];
        ry = a->ptr.p_double[p0+1];
        x0 = rx;
        y0 = ry;
        for(q=0; q<=n-2; q++)
        {
            ax = a->ptr.p_double[p0+2*kq+0];
            ay = a->ptr.p_double[p0+2*kq+1];
            buf->ptr.p_double[p1+0] = ax;
            buf->ptr.p_double[p1+1] = ay;
            rx = rx+ax;
            ry = ry+ay;
            kq = kq*rq%n;
            p1 = p1+2;
        }
        p0 = abase+aoffset+opidx*n*2;
        p1 = aoffset+opidx*n*2;
        for(q=0; q<=n-2; q++)
        {
            a->ptr.p_double[p0] = buf->ptr.p_double[p1];
            a->ptr.p_double[p0+1] = buf->ptr.p_double[p1+1];
            p0 = p0+2;
            p1 = p1+2;
        }
        
        /*
         * Convolution
         */
        ftbase_ftapplysubplan(plan, subplan, a, abase, aoffset+opidx*n*2, buf, 1, _state);
        p0 = abase+aoffset+opidx*n*2;
        p1 = precoffs;
        for(i=0; i<=n-2; i++)
        {
            ax = a->ptr.p_double[p0+0];
            ay = a->ptr.p_double[p0+1];
            bx = plan->precr.ptr.p_double[p1+0];
            by = plan->precr.ptr.p_double[p1+1];
            a->ptr.p_double[p0+0] = ax*bx-ay*by;
            a->ptr.p_double[p0+1] = -(ax*by+ay*bx);
            p0 = p0+2;
            p1 = p1+2;
        }
        ftbase_ftapplysubplan(plan, subplan, a, abase, aoffset+opidx*n*2, buf, 1, _state);
        p0 = abase+aoffset+opidx*n*2;
        for(i=0; i<=n-2; i++)
        {
            a->ptr.p_double[p0+0] = a->ptr.p_double[p0+0]/(n-1);
            a->ptr.p_double[p0+1] = -a->ptr.p_double[p0+1]/(n-1);
            p0 = p0+2;
        }
        
        /*
         * Result
         */
        buf->ptr.p_double[aoffset+opidx*n*2+0] = rx;
        buf->ptr.p_double[aoffset+opidx*n*2+1] = ry;
        kiq = 1;
        p0 = aoffset+opidx*n*2;
        p1 = abase+aoffset+opidx*n*2;
        for(q=0; q<=n-2; q++)
        {
            buf->ptr.p_double[p0+2*kiq+0] = x0+a->ptr.p_double[p1+0];
            buf->ptr.p_double[p0+2*kiq+1] = y0+a->ptr.p_double[p1+1];
            kiq = kiq*riq%n;
            p1 = p1+2;
        }
        p0 = abase+aoffset+opidx*n*2;
        p1 = aoffset+opidx*n*2;
        for(q=0; q<=n-1; q++)
        {
            a->ptr.p_double[p0] = buf->ptr.p_double[p1];
            a->ptr.p_double[p0+1] = buf->ptr.p_double[p1+1];
            p0 = p0+2;
            p1 = p1+2;
        }
    }
}


/*************************************************************************
Factorizes task size N into product of two smaller sizes N1 and N2

INPUT PARAMETERS:
    N       -   task size, N>0
    IsRoot  -   whether taks is root task (first one in a sequence)
    
OUTPUT PARAMETERS:
    N1, N2  -   such numbers that:
                * for prime N:                  N1=N2=0
                * for composite N<=MaxRadix:    N1=N2=0
                * for composite N>MaxRadix:     1<=N1<=N2, N1*N2=N

  -- ALGLIB --
     Copyright 08.04.2013 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftfactorize(ae_int_t n,
     ae_bool isroot,
     ae_int_t* n1,
     ae_int_t* n2,
     ae_state *_state)
{
    ae_int_t j;
    ae_int_t k;

    *n1 = 0;
    *n2 = 0;

    ae_assert(n>0, "FTFactorize: N<=0", _state);
    *n1 = 0;
    *n2 = 0;
    
    /*
     * Small N
     */
    if( n<=ftbase_maxradix )
    {
        return;
    }
    
    /*
     * Large N, recursive split
     */
    if( n>ftbase_recursivethreshold )
    {
        k = ae_iceil(ae_sqrt((double)(n), _state), _state)+1;
        ae_assert(k*k>=n, "FTFactorize: internal error during recursive factorization", _state);
        for(j=k; j>=2; j--)
        {
            if( n%j==0 )
            {
                *n1 = ae_minint(n/j, j, _state);
                *n2 = ae_maxint(n/j, j, _state);
                return;
            }
        }
    }
    
    /*
     * N>MaxRadix, try to find good codelet
     */
    for(j=ftbase_maxradix; j>=2; j--)
    {
        if( n%j==0 )
        {
            *n1 = j;
            *n2 = n/j;
            break;
        }
    }
    
    /*
     * In case no good codelet was found,
     * try to factorize N into product of ANY primes.
     */
    if( *n1*(*n2)!=n )
    {
        for(j=2; j<=n-1; j++)
        {
            if( n%j==0 )
            {
                *n1 = j;
                *n2 = n/j;
                break;
            }
            if( j*j>n )
            {
                break;
            }
        }
    }
    
    /*
     * normalize
     */
    if( *n1>(*n2) )
    {
        j = *n1;
        *n1 = *n2;
        *n2 = j;
    }
}


/*************************************************************************
Returns optimistic estimate of the FFT cost, in UNITs (1 UNIT = 100 KFLOPs)

INPUT PARAMETERS:
    N       -   task size, N>0
    
RESULU:
    cost in UNITs, rounded down to nearest integer

NOTE: If FFT cost is less than 1 UNIT, it will return 0 as result.

  -- ALGLIB --
     Copyright 08.04.2013 by Bochkanov Sergey
*************************************************************************/
static ae_int_t ftbase_ftoptimisticestimate(ae_int_t n, ae_state *_state)
{
    ae_int_t result;


    ae_assert(n>0, "FTOptimisticEstimate: N<=0", _state);
    result = ae_ifloor(1.0E-5*5*n*ae_log((double)(n), _state)/ae_log((double)(2), _state), _state);
    return result;
}


/*************************************************************************
Twiddle factors calculation

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ffttwcalc(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ae_int_t n1,
     ae_int_t n2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j2;
    ae_int_t n;
    ae_int_t halfn1;
    ae_int_t offs;
    double x;
    double y;
    double twxm1;
    double twy;
    double twbasexm1;
    double twbasey;
    double twrowxm1;
    double twrowy;
    double tmpx;
    double tmpy;
    double v;
    ae_int_t updatetw2;


    
    /*
     * Multiplication by twiddle factors for complex Cooley-Tukey FFT
     * with N factorized as N1*N2.
     *
     * Naive solution to this problem is given below:
     *
     *     > for K:=1 to N2-1 do
     *     >     for J:=1 to N1-1 do
     *     >     begin
     *     >         Idx:=K*N1+J;
     *     >         X:=A[AOffset+2*Idx+0];
     *     >         Y:=A[AOffset+2*Idx+1];
     *     >         TwX:=Cos(-2*Pi()*K*J/(N1*N2));
     *     >         TwY:=Sin(-2*Pi()*K*J/(N1*N2));
     *     >         A[AOffset+2*Idx+0]:=X*TwX-Y*TwY;
     *     >         A[AOffset+2*Idx+1]:=X*TwY+Y*TwX;
     *     >     end;
     *
     * However, there are exist more efficient solutions.
     *
     * Each pass of the inner cycle corresponds to multiplication of one
     * entry of A by W[k,j]=exp(-I*2*pi*k*j/N). This factor can be rewritten
     * as exp(-I*2*pi*k/N)^j. So we can replace costly exponentiation by
     * repeated multiplication: W[k,j+1]=W[k,j]*exp(-I*2*pi*k/N), with
     * second factor being computed once in the beginning of the iteration.
     *
     * Also, exp(-I*2*pi*k/N) can be represented as exp(-I*2*pi/N)^k, i.e.
     * we have W[K+1,1]=W[K,1]*W[1,1].
     *
     * In our loop we use following variables:
     * * [TwBaseXM1,TwBaseY] =   [cos(2*pi/N)-1,     sin(2*pi/N)]
     * * [TwRowXM1, TwRowY]  =   [cos(2*pi*I/N)-1,   sin(2*pi*I/N)]
     * * [TwXM1,    TwY]     =   [cos(2*pi*I*J/N)-1, sin(2*pi*I*J/N)]
     *
     * Meaning of the variables:
     * * [TwXM1,TwY] is current twiddle factor W[I,J]
     * * [TwRowXM1, TwRowY] is W[I,1]
     * * [TwBaseXM1,TwBaseY] is W[1,1]
     *
     * During inner loop we multiply current twiddle factor by W[I,1],
     * during outer loop we update W[I,1].
     *
     */
    ae_assert(ftbase_updatetw>=2, "FFTTwCalc: internal error - UpdateTw<2", _state);
    updatetw2 = ftbase_updatetw/2;
    halfn1 = n1/2;
    n = n1*n2;
    v = -2*ae_pi/n;
    twbasexm1 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
    twbasey = ae_sin(v, _state);
    twrowxm1 = (double)(0);
    twrowy = (double)(0);
    offs = aoffset;
    for(i=0; i<=n2-1; i++)
    {
        
        /*
         * Initialize twiddle factor for current row
         */
        twxm1 = (double)(0);
        twy = (double)(0);
        
        /*
         * N1-point block is separated into 2-point chunks and residual 1-point chunk
         * (in case N1 is odd). Unrolled loop is several times faster.
         */
        for(j2=0; j2<=halfn1-1; j2++)
        {
            
            /*
             * Processing:
             * * process first element in a chunk.
             * * update twiddle factor (unconditional update)
             * * process second element
             * * conditional update of the twiddle factor
             */
            x = a->ptr.p_double[offs+0];
            y = a->ptr.p_double[offs+1];
            tmpx = x*(1+twxm1)-y*twy;
            tmpy = x*twy+y*(1+twxm1);
            a->ptr.p_double[offs+0] = tmpx;
            a->ptr.p_double[offs+1] = tmpy;
            tmpx = (1+twxm1)*twrowxm1-twy*twrowy;
            twy = twy+(1+twxm1)*twrowy+twy*twrowxm1;
            twxm1 = twxm1+tmpx;
            x = a->ptr.p_double[offs+2];
            y = a->ptr.p_double[offs+3];
            tmpx = x*(1+twxm1)-y*twy;
            tmpy = x*twy+y*(1+twxm1);
            a->ptr.p_double[offs+2] = tmpx;
            a->ptr.p_double[offs+3] = tmpy;
            offs = offs+4;
            if( (j2+1)%updatetw2==0&&j2<halfn1-1 )
            {
                
                /*
                 * Recalculate twiddle factor
                 */
                v = -2*ae_pi*i*2*(j2+1)/n;
                twxm1 = ae_sin(0.5*v, _state);
                twxm1 = -2*twxm1*twxm1;
                twy = ae_sin(v, _state);
            }
            else
            {
                
                /*
                 * Update twiddle factor
                 */
                tmpx = (1+twxm1)*twrowxm1-twy*twrowy;
                twy = twy+(1+twxm1)*twrowy+twy*twrowxm1;
                twxm1 = twxm1+tmpx;
            }
        }
        if( n1%2==1 )
        {
            
            /*
             * Handle residual chunk
             */
            x = a->ptr.p_double[offs+0];
            y = a->ptr.p_double[offs+1];
            tmpx = x*(1+twxm1)-y*twy;
            tmpy = x*twy+y*(1+twxm1);
            a->ptr.p_double[offs+0] = tmpx;
            a->ptr.p_double[offs+1] = tmpy;
            offs = offs+2;
        }
        
        /*
         * update TwRow: TwRow(new) = TwRow(old)*TwBase
         */
        if( i<n2-1 )
        {
            if( (i+1)%ftbase_updatetw==0 )
            {
                v = -2*ae_pi*(i+1)/n;
                twrowxm1 = ae_sin(0.5*v, _state);
                twrowxm1 = -2*twrowxm1*twrowxm1;
                twrowy = ae_sin(v, _state);
            }
            else
            {
                tmpx = twbasexm1+twrowxm1*twbasexm1-twrowy*twbasey;
                tmpy = twbasey+twrowxm1*twbasey+twrowy*twbasexm1;
                twrowxm1 = twrowxm1+tmpx;
                twrowy = twrowy+tmpy;
            }
        }
    }
}


/*************************************************************************
Linear transpose: transpose complex matrix stored in 1-dimensional array

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_internalcomplexlintranspose(/* Real    */ ae_vector* a,
     ae_int_t m,
     ae_int_t n,
     ae_int_t astart,
     /* Real    */ ae_vector* buf,
     ae_state *_state)
{


    ftbase_ffticltrec(a, astart, n, buf, 0, m, m, n, _state);
    ae_v_move(&a->ptr.p_double[astart], 1, &buf->ptr.p_double[0], 1, ae_v_len(astart,astart+2*m*n-1));
}


/*************************************************************************
Recurrent subroutine for a InternalComplexLinTranspose

Write A^T to B, where:
* A is m*n complex matrix stored in array A as pairs of real/image values,
  beginning from AStart position, with AStride stride
* B is n*m complex matrix stored in array B as pairs of real/image values,
  beginning from BStart position, with BStride stride
stride is measured in complex numbers, i.e. in real/image pairs.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ffticltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t idx1;
    ae_int_t idx2;
    ae_int_t m2;
    ae_int_t m1;
    ae_int_t n1;


    if( m==0||n==0 )
    {
        return;
    }
    if( ae_maxint(m, n, _state)<=8 )
    {
        m2 = 2*bstride;
        for(i=0; i<=m-1; i++)
        {
            idx1 = bstart+2*i;
            idx2 = astart+2*i*astride;
            for(j=0; j<=n-1; j++)
            {
                b->ptr.p_double[idx1+0] = a->ptr.p_double[idx2+0];
                b->ptr.p_double[idx1+1] = a->ptr.p_double[idx2+1];
                idx1 = idx1+m2;
                idx2 = idx2+2;
            }
        }
        return;
    }
    if( n>m )
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
         *                                  ( B2 )
         */
        n1 = n/2;
        if( n-n1>=8&&n1%8!=0 )
        {
            n1 = n1+(8-n1%8);
        }
        ae_assert(n-n1>0, "Assertion failed", _state);
        ftbase_ffticltrec(a, astart, astride, b, bstart, bstride, m, n1, _state);
        ftbase_ffticltrec(a, astart+2*n1, astride, b, bstart+2*n1*bstride, bstride, m, n-n1, _state);
    }
    else
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
         *                     ( A2 )
         */
        m1 = m/2;
        if( m-m1>=8&&m1%8!=0 )
        {
            m1 = m1+(8-m1%8);
        }
        ae_assert(m-m1>0, "Assertion failed", _state);
        ftbase_ffticltrec(a, astart, astride, b, bstart, bstride, m1, n, _state);
        ftbase_ffticltrec(a, astart+2*m1*astride, astride, b, bstart+2*m1, bstride, m-m1, n, _state);
    }
}


/*************************************************************************
Recurrent subroutine for a InternalRealLinTranspose


  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_fftirltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t idx1;
    ae_int_t idx2;
    ae_int_t m1;
    ae_int_t n1;


    if( m==0||n==0 )
    {
        return;
    }
    if( ae_maxint(m, n, _state)<=8 )
    {
        for(i=0; i<=m-1; i++)
        {
            idx1 = bstart+i;
            idx2 = astart+i*astride;
            for(j=0; j<=n-1; j++)
            {
                b->ptr.p_double[idx1] = a->ptr.p_double[idx2];
                idx1 = idx1+bstride;
                idx2 = idx2+1;
            }
        }
        return;
    }
    if( n>m )
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
         *                                  ( B2 )
         */
        n1 = n/2;
        if( n-n1>=8&&n1%8!=0 )
        {
            n1 = n1+(8-n1%8);
        }
        ae_assert(n-n1>0, "Assertion failed", _state);
        ftbase_fftirltrec(a, astart, astride, b, bstart, bstride, m, n1, _state);
        ftbase_fftirltrec(a, astart+n1, astride, b, bstart+n1*bstride, bstride, m, n-n1, _state);
    }
    else
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
         *                     ( A2 )
         */
        m1 = m/2;
        if( m-m1>=8&&m1%8!=0 )
        {
            m1 = m1+(8-m1%8);
        }
        ae_assert(m-m1>0, "Assertion failed", _state);
        ftbase_fftirltrec(a, astart, astride, b, bstart, bstride, m1, n, _state);
        ftbase_fftirltrec(a, astart+m1*astride, astride, b, bstart+m1, bstride, m-m1, n, _state);
    }
}


/*************************************************************************
recurrent subroutine for FFTFindSmoothRec

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftbasefindsmoothrec(ae_int_t n,
     ae_int_t seed,
     ae_int_t leastfactor,
     ae_int_t* best,
     ae_state *_state)
{


    ae_assert(ftbase_ftbasemaxsmoothfactor<=5, "FTBaseFindSmoothRec: internal error!", _state);
    if( seed>=n )
    {
        *best = ae_minint(*best, seed, _state);
        return;
    }
    if( leastfactor<=2 )
    {
        ftbase_ftbasefindsmoothrec(n, seed*2, 2, best, _state);
    }
    if( leastfactor<=3 )
    {
        ftbase_ftbasefindsmoothrec(n, seed*3, 3, best, _state);
    }
    if( leastfactor<=5 )
    {
        ftbase_ftbasefindsmoothrec(n, seed*5, 5, best, _state);
    }
}


void _fasttransformplan_init(void* _p, ae_state *_state)
{
    fasttransformplan *p = (fasttransformplan*)_p;
    ae_touch_ptr((void*)p);
    ae_matrix_init(&p->entries, 0, 0, DT_INT, _state);
    ae_vector_init(&p->buffer, 0, DT_REAL, _state);
    ae_vector_init(&p->precr, 0, DT_REAL, _state);
    ae_vector_init(&p->preci, 0, DT_REAL, _state);
    ae_shared_pool_init(&p->bluesteinpool, _state);
}


void _fasttransformplan_init_copy(void* _dst, void* _src, ae_state *_state)
{
    fasttransformplan *dst = (fasttransformplan*)_dst;
    fasttransformplan *src = (fasttransformplan*)_src;
    ae_matrix_init_copy(&dst->entries, &src->entries, _state);
    ae_vector_init_copy(&dst->buffer, &src->buffer, _state);
    ae_vector_init_copy(&dst->precr, &src->precr, _state);
    ae_vector_init_copy(&dst->preci, &src->preci, _state);
    ae_shared_pool_init_copy(&dst->bluesteinpool, &src->bluesteinpool, _state);
}


void _fasttransformplan_clear(void* _p)
{
    fasttransformplan *p = (fasttransformplan*)_p;
    ae_touch_ptr((void*)p);
    ae_matrix_clear(&p->entries);
    ae_vector_clear(&p->buffer);
    ae_vector_clear(&p->precr);
    ae_vector_clear(&p->preci);
    ae_shared_pool_clear(&p->bluesteinpool);
}


void _fasttransformplan_destroy(void* _p)
{
    fasttransformplan *p = (fasttransformplan*)_p;
    ae_touch_ptr((void*)p);
    ae_matrix_destroy(&p->entries);
    ae_vector_destroy(&p->buffer);
    ae_vector_destroy(&p->precr);
    ae_vector_destroy(&p->preci);
    ae_shared_pool_destroy(&p->bluesteinpool);
}




double nulog1p(double x, ae_state *_state)
{
    double z;
    double lp;
    double lq;
    double result;


    z = 1.0+x;
    if( ae_fp_less(z,0.70710678118654752440)||ae_fp_greater(z,1.41421356237309504880) )
    {
        result = ae_log(z, _state);
        return result;
    }
    z = x*x;
    lp = 4.5270000862445199635215E-5;
    lp = lp*x+4.9854102823193375972212E-1;
    lp = lp*x+6.5787325942061044846969E0;
    lp = lp*x+2.9911919328553073277375E1;
    lp = lp*x+6.0949667980987787057556E1;
    lp = lp*x+5.7112963590585538103336E1;
    lp = lp*x+2.0039553499201281259648E1;
    lq = 1.0000000000000000000000E0;
    lq = lq*x+1.5062909083469192043167E1;
    lq = lq*x+8.3047565967967209469434E1;
    lq = lq*x+2.2176239823732856465394E2;
    lq = lq*x+3.0909872225312059774938E2;
    lq = lq*x+2.1642788614495947685003E2;
    lq = lq*x+6.0118660497603843919306E1;
    z = -0.5*z+x*(z*lp/lq);
    result = x+z;
    return result;
}


double nuexpm1(double x, ae_state *_state)
{
    double r;
    double xx;
    double ep;
    double eq;
    double result;


    if( ae_fp_less(x,-0.5)||ae_fp_greater(x,0.5) )
    {
        result = ae_exp(x, _state)-1.0;
        return result;
    }
    xx = x*x;
    ep = 1.2617719307481059087798E-4;
    ep = ep*xx+3.0299440770744196129956E-2;
    ep = ep*xx+9.9999999999999999991025E-1;
    eq = 3.0019850513866445504159E-6;
    eq = eq*xx+2.5244834034968410419224E-3;
    eq = eq*xx+2.2726554820815502876593E-1;
    eq = eq*xx+2.0000000000000000000897E0;
    r = x*ep;
    r = r/(eq-r);
    result = r+r;
    return result;
}


double nucosm1(double x, ae_state *_state)
{
    double xx;
    double c;
    double result;


    if( ae_fp_less(x,-0.25*ae_pi)||ae_fp_greater(x,0.25*ae_pi) )
    {
        result = ae_cos(x, _state)-1;
        return result;
    }
    xx = x*x;
    c = 4.7377507964246204691685E-14;
    c = c*xx-1.1470284843425359765671E-11;
    c = c*xx+2.0876754287081521758361E-9;
    c = c*xx-2.7557319214999787979814E-7;
    c = c*xx+2.4801587301570552304991E-5;
    c = c*xx-1.3888888888888872993737E-3;
    c = c*xx+4.1666666666666666609054E-2;
    result = -0.5*xx+xx*xx*c;
    return result;
}





}

