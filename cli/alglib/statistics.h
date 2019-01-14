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
#ifndef _statistics_pkg_h
#define _statistics_pkg_h
#include "ap.h"
#include "alglibinternal.h"
#include "linalg.h"
#include "specialfunctions.h"

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
Calculation of the distribution moments: mean, variance, skewness, kurtosis.

INPUT PARAMETERS:
    X       -   sample
    N       -   N>=0, sample size:
                * if given, only leading N elements of X are processed
                * if not given, automatically determined from size of X

OUTPUT PARAMETERS
    Mean    -   mean.
    Variance-   variance.
    Skewness-   skewness (if variance<>0; zero otherwise).
    Kurtosis-   kurtosis (if variance<>0; zero otherwise).

NOTE: variance is calculated by dividing sum of squares by N-1, not N.

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
void samplemoments(const real_1d_array &x, const ae_int_t n, double &mean, double &variance, double &skewness, double &kurtosis);
void samplemoments(const real_1d_array &x, double &mean, double &variance, double &skewness, double &kurtosis);


/*************************************************************************
Calculation of the mean.

INPUT PARAMETERS:
    X       -   sample
    N       -   N>=0, sample size:
                * if given, only leading N elements of X are processed
                * if not given, automatically determined from size of X

NOTE:

This function return result  which calculated by 'SampleMoments' function
and stored at 'Mean' variable.


  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
double samplemean(const real_1d_array &x, const ae_int_t n);
double samplemean(const real_1d_array &x);


/*************************************************************************
Calculation of the variance.

INPUT PARAMETERS:
    X       -   sample
    N       -   N>=0, sample size:
                * if given, only leading N elements of X are processed
                * if not given, automatically determined from size of X

NOTE:

This function return result  which calculated by 'SampleMoments' function
and stored at 'Variance' variable.


  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
double samplevariance(const real_1d_array &x, const ae_int_t n);
double samplevariance(const real_1d_array &x);


/*************************************************************************
Calculation of the skewness.

INPUT PARAMETERS:
    X       -   sample
    N       -   N>=0, sample size:
                * if given, only leading N elements of X are processed
                * if not given, automatically determined from size of X

NOTE:

This function return result  which calculated by 'SampleMoments' function
and stored at 'Skewness' variable.


  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
double sampleskewness(const real_1d_array &x, const ae_int_t n);
double sampleskewness(const real_1d_array &x);


/*************************************************************************
Calculation of the kurtosis.

INPUT PARAMETERS:
    X       -   sample
    N       -   N>=0, sample size:
                * if given, only leading N elements of X are processed
                * if not given, automatically determined from size of X

NOTE:

This function return result  which calculated by 'SampleMoments' function
and stored at 'Kurtosis' variable.


  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
double samplekurtosis(const real_1d_array &x, const ae_int_t n);
double samplekurtosis(const real_1d_array &x);


/*************************************************************************
ADev

Input parameters:
    X   -   sample
    N   -   N>=0, sample size:
            * if given, only leading N elements of X are processed
            * if not given, automatically determined from size of X

Output parameters:
    ADev-   ADev

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
void sampleadev(const real_1d_array &x, const ae_int_t n, double &adev);
void sampleadev(const real_1d_array &x, double &adev);


/*************************************************************************
Median calculation.

Input parameters:
    X   -   sample (array indexes: [0..N-1])
    N   -   N>=0, sample size:
            * if given, only leading N elements of X are processed
            * if not given, automatically determined from size of X

Output parameters:
    Median

  -- ALGLIB --
     Copyright 06.09.2006 by Bochkanov Sergey
*************************************************************************/
void samplemedian(const real_1d_array &x, const ae_int_t n, double &median);
void samplemedian(const real_1d_array &x, double &median);


/*************************************************************************
Percentile calculation.

Input parameters:
    X   -   sample (array indexes: [0..N-1])
    N   -   N>=0, sample size:
            * if given, only leading N elements of X are processed
            * if not given, automatically determined from size of X
    P   -   percentile (0<=P<=1)

Output parameters:
    V   -   percentile

  -- ALGLIB --
     Copyright 01.03.2008 by Bochkanov Sergey
*************************************************************************/
void samplepercentile(const real_1d_array &x, const ae_int_t n, const double p, double &v);
void samplepercentile(const real_1d_array &x, const double p, double &v);


/*************************************************************************
2-sample covariance

Input parameters:
    X       -   sample 1 (array indexes: [0..N-1])
    Y       -   sample 2 (array indexes: [0..N-1])
    N       -   N>=0, sample size:
                * if given, only N leading elements of X/Y are processed
                * if not given, automatically determined from input sizes

Result:
    covariance (zero for N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
double cov2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
double cov2(const real_1d_array &x, const real_1d_array &y);


/*************************************************************************
Pearson product-moment correlation coefficient

Input parameters:
    X       -   sample 1 (array indexes: [0..N-1])
    Y       -   sample 2 (array indexes: [0..N-1])
    N       -   N>=0, sample size:
                * if given, only N leading elements of X/Y are processed
                * if not given, automatically determined from input sizes

Result:
    Pearson product-moment correlation coefficient
    (zero for N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
double pearsoncorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
double pearsoncorr2(const real_1d_array &x, const real_1d_array &y);


/*************************************************************************
Spearman's rank correlation coefficient

Input parameters:
    X       -   sample 1 (array indexes: [0..N-1])
    Y       -   sample 2 (array indexes: [0..N-1])
    N       -   N>=0, sample size:
                * if given, only N leading elements of X/Y are processed
                * if not given, automatically determined from input sizes

Result:
    Spearman's rank correlation coefficient
    (zero for N=0 or N=1)

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
double spearmancorr2(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);
double spearmancorr2(const real_1d_array &x, const real_1d_array &y);


/*************************************************************************
Covariance matrix

SMP EDITION OF ALGLIB:

  ! This function can utilize multicore capabilities of  your system.  In
  ! order to do this you have to call version with "smp_" prefix,   which
  ! indicates that multicore code will be used.
  !
  ! This note is given for users of SMP edition; if you use GPL  edition,
  ! or commercial edition of ALGLIB without SMP support, you  still  will
  ! be able to call smp-version of this function,  but  all  computations
  ! will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.
  !
  ! You should remember that starting/stopping worker thread always  have
  ! non-zero cost. Although  multicore  version  is  pretty  efficient on
  ! large problems, we do not recommend you to use it on small problems -
  ! with covariance matrices smaller than 128*128.

INPUT PARAMETERS:
    X   -   array[N,M], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X are used
            * if not given, automatically determined from input size
    M   -   M>0, number of variables:
            * if given, only leading M columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M,M], covariance matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void covm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
void smp_covm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
void covm(const real_2d_array &x, real_2d_array &c);
void smp_covm(const real_2d_array &x, real_2d_array &c);


/*************************************************************************
Pearson product-moment correlation matrix

SMP EDITION OF ALGLIB:

  ! This function can utilize multicore capabilities of  your system.  In
  ! order to do this you have to call version with "smp_" prefix,   which
  ! indicates that multicore code will be used.
  !
  ! This note is given for users of SMP edition; if you use GPL  edition,
  ! or commercial edition of ALGLIB without SMP support, you  still  will
  ! be able to call smp-version of this function,  but  all  computations
  ! will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.
  !
  ! You should remember that starting/stopping worker thread always  have
  ! non-zero cost. Although  multicore  version  is  pretty  efficient on
  ! large problems, we do not recommend you to use it on small problems -
  ! with correlation matrices smaller than 128*128.

INPUT PARAMETERS:
    X   -   array[N,M], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X are used
            * if not given, automatically determined from input size
    M   -   M>0, number of variables:
            * if given, only leading M columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M,M], correlation matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void pearsoncorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
void smp_pearsoncorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
void pearsoncorrm(const real_2d_array &x, real_2d_array &c);
void smp_pearsoncorrm(const real_2d_array &x, real_2d_array &c);


/*************************************************************************
Spearman's rank correlation matrix

SMP EDITION OF ALGLIB:

  ! This function can utilize multicore capabilities of  your system.  In
  ! order to do this you have to call version with "smp_" prefix,   which
  ! indicates that multicore code will be used.
  !
  ! This note is given for users of SMP edition; if you use GPL  edition,
  ! or commercial edition of ALGLIB without SMP support, you  still  will
  ! be able to call smp-version of this function,  but  all  computations
  ! will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.
  !
  ! You should remember that starting/stopping worker thread always  have
  ! non-zero cost. Although  multicore  version  is  pretty  efficient on
  ! large problems, we do not recommend you to use it on small problems -
  ! with correlation matrices smaller than 128*128.

INPUT PARAMETERS:
    X   -   array[N,M], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X are used
            * if not given, automatically determined from input size
    M   -   M>0, number of variables:
            * if given, only leading M columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M,M], correlation matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void spearmancorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
void smp_spearmancorrm(const real_2d_array &x, const ae_int_t n, const ae_int_t m, real_2d_array &c);
void spearmancorrm(const real_2d_array &x, real_2d_array &c);
void smp_spearmancorrm(const real_2d_array &x, real_2d_array &c);


/*************************************************************************
Cross-covariance matrix

SMP EDITION OF ALGLIB:

  ! This function can utilize multicore capabilities of  your system.  In
  ! order to do this you have to call version with "smp_" prefix,   which
  ! indicates that multicore code will be used.
  !
  ! This note is given for users of SMP edition; if you use GPL  edition,
  ! or commercial edition of ALGLIB without SMP support, you  still  will
  ! be able to call smp-version of this function,  but  all  computations
  ! will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.
  !
  ! You should remember that starting/stopping worker thread always  have
  ! non-zero cost. Although  multicore  version  is  pretty  efficient on
  ! large problems, we do not recommend you to use it on small problems -
  ! with covariance matrices smaller than 128*128.

INPUT PARAMETERS:
    X   -   array[N,M1], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    Y   -   array[N,M2], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X/Y are used
            * if not given, automatically determined from input sizes
    M1  -   M1>0, number of variables in X:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size
    M2  -   M2>0, number of variables in Y:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M1,M2], cross-covariance matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void covm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
void smp_covm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
void covm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
void smp_covm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);


/*************************************************************************
Pearson product-moment cross-correlation matrix

SMP EDITION OF ALGLIB:

  ! This function can utilize multicore capabilities of  your system.  In
  ! order to do this you have to call version with "smp_" prefix,   which
  ! indicates that multicore code will be used.
  !
  ! This note is given for users of SMP edition; if you use GPL  edition,
  ! or commercial edition of ALGLIB without SMP support, you  still  will
  ! be able to call smp-version of this function,  but  all  computations
  ! will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.
  !
  ! You should remember that starting/stopping worker thread always  have
  ! non-zero cost. Although  multicore  version  is  pretty  efficient on
  ! large problems, we do not recommend you to use it on small problems -
  ! with correlation matrices smaller than 128*128.

INPUT PARAMETERS:
    X   -   array[N,M1], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    Y   -   array[N,M2], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X/Y are used
            * if not given, automatically determined from input sizes
    M1  -   M1>0, number of variables in X:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size
    M2  -   M2>0, number of variables in Y:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M1,M2], cross-correlation matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
void smp_pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
void pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
void smp_pearsoncorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);


/*************************************************************************
Spearman's rank cross-correlation matrix

SMP EDITION OF ALGLIB:

  ! This function can utilize multicore capabilities of  your system.  In
  ! order to do this you have to call version with "smp_" prefix,   which
  ! indicates that multicore code will be used.
  !
  ! This note is given for users of SMP edition; if you use GPL  edition,
  ! or commercial edition of ALGLIB without SMP support, you  still  will
  ! be able to call smp-version of this function,  but  all  computations
  ! will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.
  !
  ! You should remember that starting/stopping worker thread always  have
  ! non-zero cost. Although  multicore  version  is  pretty  efficient on
  ! large problems, we do not recommend you to use it on small problems -
  ! with correlation matrices smaller than 128*128.

INPUT PARAMETERS:
    X   -   array[N,M1], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    Y   -   array[N,M2], sample matrix:
            * J-th column corresponds to J-th variable
            * I-th row corresponds to I-th observation
    N   -   N>=0, number of observations:
            * if given, only leading N rows of X/Y are used
            * if not given, automatically determined from input sizes
    M1  -   M1>0, number of variables in X:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size
    M2  -   M2>0, number of variables in Y:
            * if given, only leading M1 columns of X are used
            * if not given, automatically determined from input size

OUTPUT PARAMETERS:
    C   -   array[M1,M2], cross-correlation matrix (zero if N=0 or N=1)

  -- ALGLIB --
     Copyright 28.10.2010 by Bochkanov Sergey
*************************************************************************/
void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
void smp_spearmancorrm2(const real_2d_array &x, const real_2d_array &y, const ae_int_t n, const ae_int_t m1, const ae_int_t m2, real_2d_array &c);
void spearmancorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);
void smp_spearmancorrm2(const real_2d_array &x, const real_2d_array &y, real_2d_array &c);


/*************************************************************************
This function replaces data in XY by their ranks:
* XY is processed row-by-row
* rows are processed separately
* tied data are correctly handled (tied ranks are calculated)
* ranking starts from 0, ends at NFeatures-1
* sum of within-row values is equal to (NFeatures-1)*NFeatures/2

SMP EDITION OF ALGLIB:

  ! This function can utilize multicore capabilities of  your system.  In
  ! order to do this you have to call version with "smp_" prefix,   which
  ! indicates that multicore code will be used.
  !
  ! This note is given for users of SMP edition; if you use GPL  edition,
  ! or commercial edition of ALGLIB without SMP support, you  still  will
  ! be able to call smp-version of this function,  but  all  computations
  ! will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.
  !
  ! You should remember that starting/stopping worker thread always  have
  ! non-zero cost. Although  multicore  version  is  pretty  efficient on
  ! large problems, we do not recommend you to use it on small problems -
  ! ones where expected operations count is less than 100.000

INPUT PARAMETERS:
    XY      -   array[NPoints,NFeatures], dataset
    NPoints -   number of points
    NFeatures-  number of features

OUTPUT PARAMETERS:
    XY      -   data are replaced by their within-row ranks;
                ranking starts from 0, ends at NFeatures-1

  -- ALGLIB --
     Copyright 18.04.2013 by Bochkanov Sergey
*************************************************************************/
void rankdata(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
void smp_rankdata(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
void rankdata(real_2d_array &xy);
void smp_rankdata(real_2d_array &xy);


/*************************************************************************
This function replaces data in XY by their CENTERED ranks:
* XY is processed row-by-row
* rows are processed separately
* tied data are correctly handled (tied ranks are calculated)
* centered ranks are just usual ranks, but centered in such way  that  sum
  of within-row values is equal to 0.0.
* centering is performed by subtracting mean from each row, i.e it changes
  mean value, but does NOT change higher moments

SMP EDITION OF ALGLIB:

  ! This function can utilize multicore capabilities of  your system.  In
  ! order to do this you have to call version with "smp_" prefix,   which
  ! indicates that multicore code will be used.
  !
  ! This note is given for users of SMP edition; if you use GPL  edition,
  ! or commercial edition of ALGLIB without SMP support, you  still  will
  ! be able to call smp-version of this function,  but  all  computations
  ! will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.
  !
  ! You should remember that starting/stopping worker thread always  have
  ! non-zero cost. Although  multicore  version  is  pretty  efficient on
  ! large problems, we do not recommend you to use it on small problems -
  ! ones where expected operations count is less than 100.000

INPUT PARAMETERS:
    XY      -   array[NPoints,NFeatures], dataset
    NPoints -   number of points
    NFeatures-  number of features

OUTPUT PARAMETERS:
    XY      -   data are replaced by their within-row ranks;
                ranking starts from 0, ends at NFeatures-1

  -- ALGLIB --
     Copyright 18.04.2013 by Bochkanov Sergey
*************************************************************************/
void rankdatacentered(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
void smp_rankdatacentered(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures);
void rankdatacentered(real_2d_array &xy);
void smp_rankdatacentered(real_2d_array &xy);


/*************************************************************************
Obsolete function, we recommend to use PearsonCorr2().

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
double pearsoncorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);


/*************************************************************************
Obsolete function, we recommend to use SpearmanCorr2().

    -- ALGLIB --
    Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
double spearmanrankcorrelation(const real_1d_array &x, const real_1d_array &y, const ae_int_t n);

/*************************************************************************
Pearson's correlation coefficient significance test

This test checks hypotheses about whether X  and  Y  are  samples  of  two
continuous  distributions  having  zero  correlation  or   whether   their
correlation is non-zero.

The following tests are performed:
    * two-tailed test (null hypothesis - X and Y have zero correlation)
    * left-tailed test (null hypothesis - the correlation  coefficient  is
      greater than or equal to 0)
    * right-tailed test (null hypothesis - the correlation coefficient  is
      less than or equal to 0).

Requirements:
    * the number of elements in each sample is not less than 5
    * normality of distributions of X and Y.

Input parameters:
    R   -   Pearson's correlation coefficient for X and Y
    N   -   number of elements in samples, N>=5.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
void pearsoncorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail);


/*************************************************************************
Spearman's rank correlation coefficient significance test

This test checks hypotheses about whether X  and  Y  are  samples  of  two
continuous  distributions  having  zero  correlation  or   whether   their
correlation is non-zero.

The following tests are performed:
    * two-tailed test (null hypothesis - X and Y have zero correlation)
    * left-tailed test (null hypothesis - the correlation  coefficient  is
      greater than or equal to 0)
    * right-tailed test (null hypothesis - the correlation coefficient  is
      less than or equal to 0).

Requirements:
    * the number of elements in each sample is not less than 5.

The test is non-parametric and doesn't require distributions X and Y to be
normal.

Input parameters:
    R   -   Spearman's rank correlation coefficient for X and Y
    N   -   number of elements in samples, N>=5.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
void spearmanrankcorrelationsignificance(const double r, const ae_int_t n, double &bothtails, double &lefttail, double &righttail);

/*************************************************************************
Jarque-Bera test

This test checks hypotheses about the fact that a  given  sample  X  is  a
sample of normal random variable.

Requirements:
    * the number of elements in the sample is not less than 5.

Input parameters:
    X   -   sample. Array whose index goes from 0 to N-1.
    N   -   size of the sample. N>=5

Output parameters:
    P           -   p-value for the test

Accuracy of the approximation used (5<=N<=1951):

p-value  	    relative error (5<=N<=1951)
[1, 0.1]            < 1%
[0.1, 0.01]         < 2%
[0.01, 0.001]       < 6%
[0.001, 0]          wasn't measured

For N>1951 accuracy wasn't measured but it shouldn't be sharply  different
from table values.

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
void jarqueberatest(const real_1d_array &x, const ae_int_t n, double &p);

/*************************************************************************
Mann-Whitney U-test

This test checks hypotheses about whether X  and  Y  are  samples  of  two
continuous distributions of the same shape  and  same  median  or  whether
their medians are different.

The following tests are performed:
    * two-tailed test (null hypothesis - the medians are equal)
    * left-tailed test (null hypothesis - the median of the  first  sample
      is greater than or equal to the median of the second sample)
    * right-tailed test (null hypothesis - the median of the first  sample
      is less than or equal to the median of the second sample).

Requirements:
    * the samples are independent
    * X and Y are continuous distributions (or discrete distributions well-
      approximating continuous distributions)
    * distributions of X and Y have the  same  shape.  The  only  possible
      difference is their position (i.e. the value of the median)
    * the number of elements in each sample is not less than 5
    * the scale of measurement should be ordinal, interval or ratio  (i.e.
      the test could not be applied to nominal variables).

The test is non-parametric and doesn't require distributions to be normal.

Input parameters:
    X   -   sample 1. Array whose index goes from 0 to N-1.
    N   -   size of the sample. N>=5
    Y   -   sample 2. Array whose index goes from 0 to M-1.
    M   -   size of the sample. M>=5

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

To calculate p-values, special approximation is used. This method lets  us
calculate p-values with satisfactory  accuracy  in  interval  [0.0001, 1].
There is no approximation outside the [0.0001, 1] interval. Therefore,  if
the significance level outlies this interval, the test returns 0.0001.

Relative precision of approximation of p-value:

N          M          Max.err.   Rms.err.
5..10      N..10      1.4e-02    6.0e-04
5..10      N..100     2.2e-02    5.3e-06
10..15     N..15      1.0e-02    3.2e-04
10..15     N..100     1.0e-02    2.2e-05
15..100    N..100     6.1e-03    2.7e-06

For N,M>100 accuracy checks weren't put into  practice,  but  taking  into
account characteristics of asymptotic approximation used, precision should
not be sharply different from the values for interval [5, 100].

NOTE: P-value approximation was  optimized  for  0.0001<=p<=0.2500.  Thus,
      P's outside of this interval are enforced to these bounds. Say,  you
      may quite often get P equal to exactly 0.25 or 0.0001.

  -- ALGLIB --
     Copyright 09.04.2007 by Bochkanov Sergey
*************************************************************************/
void mannwhitneyutest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);

/*************************************************************************
Sign test

This test checks three hypotheses about the median of  the  given  sample.
The following tests are performed:
    * two-tailed test (null hypothesis - the median is equal to the  given
      value)
    * left-tailed test (null hypothesis - the median is  greater  than  or
      equal to the given value)
    * right-tailed test (null hypothesis - the  median  is  less  than  or
      equal to the given value)

Requirements:
    * the scale of measurement should be ordinal, interval or ratio  (i.e.
      the test could not be applied to nominal variables).

The test is non-parametric and doesn't require distribution X to be normal

Input parameters:
    X       -   sample. Array whose index goes from 0 to N-1.
    N       -   size of the sample.
    Median  -   assumed median value.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

While   calculating   p-values   high-precision   binomial    distribution
approximation is used, so significance levels have about 15 exact digits.

  -- ALGLIB --
     Copyright 08.09.2006 by Bochkanov Sergey
*************************************************************************/
void onesamplesigntest(const real_1d_array &x, const ae_int_t n, const double median, double &bothtails, double &lefttail, double &righttail);

/*************************************************************************
One-sample t-test

This test checks three hypotheses about the mean of the given sample.  The
following tests are performed:
    * two-tailed test (null hypothesis - the mean is equal  to  the  given
      value)
    * left-tailed test (null hypothesis - the  mean  is  greater  than  or
      equal to the given value)
    * right-tailed test (null hypothesis - the mean is less than or  equal
      to the given value).

The test is based on the assumption that  a  given  sample  has  a  normal
distribution and  an  unknown  dispersion.  If  the  distribution  sharply
differs from normal, the test will work incorrectly.

INPUT PARAMETERS:
    X       -   sample. Array whose index goes from 0 to N-1.
    N       -   size of sample, N>=0
    Mean    -   assumed value of the mean.

OUTPUT PARAMETERS:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

NOTE: this function correctly handles degenerate cases:
      * when N=0, all p-values are set to 1.0
      * when variance of X[] is exactly zero, p-values are set
        to 1.0 or 0.0, depending on difference between sample mean and
        value of mean being tested.


  -- ALGLIB --
     Copyright 08.09.2006 by Bochkanov Sergey
*************************************************************************/
void studentttest1(const real_1d_array &x, const ae_int_t n, const double mean, double &bothtails, double &lefttail, double &righttail);


/*************************************************************************
Two-sample pooled test

This test checks three hypotheses about the mean of the given samples. The
following tests are performed:
    * two-tailed test (null hypothesis - the means are equal)
    * left-tailed test (null hypothesis - the mean of the first sample  is
      greater than or equal to the mean of the second sample)
    * right-tailed test (null hypothesis - the mean of the first sample is
      less than or equal to the mean of the second sample).

Test is based on the following assumptions:
    * given samples have normal distributions
    * dispersions are equal
    * samples are independent.

Input parameters:
    X       -   sample 1. Array whose index goes from 0 to N-1.
    N       -   size of sample.
    Y       -   sample 2. Array whose index goes from 0 to M-1.
    M       -   size of sample.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

NOTE: this function correctly handles degenerate cases:
      * when N=0 or M=0, all p-values are set to 1.0
      * when both samples has exactly zero variance, p-values are set
        to 1.0 or 0.0, depending on difference between means.

  -- ALGLIB --
     Copyright 18.09.2006 by Bochkanov Sergey
*************************************************************************/
void studentttest2(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);


/*************************************************************************
Two-sample unpooled test

This test checks three hypotheses about the mean of the given samples. The
following tests are performed:
    * two-tailed test (null hypothesis - the means are equal)
    * left-tailed test (null hypothesis - the mean of the first sample  is
      greater than or equal to the mean of the second sample)
    * right-tailed test (null hypothesis - the mean of the first sample is
      less than or equal to the mean of the second sample).

Test is based on the following assumptions:
    * given samples have normal distributions
    * samples are independent.
Equality of variances is NOT required.

Input parameters:
    X - sample 1. Array whose index goes from 0 to N-1.
    N - size of the sample.
    Y - sample 2. Array whose index goes from 0 to M-1.
    M - size of the sample.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

NOTE: this function correctly handles degenerate cases:
      * when N=0 or M=0, all p-values are set to 1.0
      * when both samples has zero variance, p-values are set
        to 1.0 or 0.0, depending on difference between means.
      * when only one sample has zero variance, test reduces to 1-sample
        version.

  -- ALGLIB --
     Copyright 18.09.2006 by Bochkanov Sergey
*************************************************************************/
void unequalvariancettest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);

/*************************************************************************
Two-sample F-test

This test checks three hypotheses about dispersions of the given  samples.
The following tests are performed:
    * two-tailed test (null hypothesis - the dispersions are equal)
    * left-tailed test (null hypothesis  -  the  dispersion  of  the first
      sample is greater than or equal to  the  dispersion  of  the  second
      sample).
    * right-tailed test (null hypothesis - the  dispersion  of  the  first
      sample is less than or equal to the dispersion of the second sample)

The test is based on the following assumptions:
    * the given samples have normal distributions
    * the samples are independent.

Input parameters:
    X   -   sample 1. Array whose index goes from 0 to N-1.
    N   -   sample size.
    Y   -   sample 2. Array whose index goes from 0 to M-1.
    M   -   sample size.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 19.09.2006 by Bochkanov Sergey
*************************************************************************/
void ftest(const real_1d_array &x, const ae_int_t n, const real_1d_array &y, const ae_int_t m, double &bothtails, double &lefttail, double &righttail);


/*************************************************************************
One-sample chi-square test

This test checks three hypotheses about the dispersion of the given sample
The following tests are performed:
    * two-tailed test (null hypothesis - the dispersion equals  the  given
      number)
    * left-tailed test (null hypothesis - the dispersion is  greater  than
      or equal to the given number)
    * right-tailed test (null hypothesis  -  dispersion is  less  than  or
      equal to the given number).

Test is based on the following assumptions:
    * the given sample has a normal distribution.

Input parameters:
    X           -   sample 1. Array whose index goes from 0 to N-1.
    N           -   size of the sample.
    Variance    -   dispersion value to compare with.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

  -- ALGLIB --
     Copyright 19.09.2006 by Bochkanov Sergey
*************************************************************************/
void onesamplevariancetest(const real_1d_array &x, const ae_int_t n, const double variance, double &bothtails, double &lefttail, double &righttail);

/*************************************************************************
Wilcoxon signed-rank test

This test checks three hypotheses about the median  of  the  given sample.
The following tests are performed:
    * two-tailed test (null hypothesis - the median is equal to the  given
      value)
    * left-tailed test (null hypothesis - the median is  greater  than  or
      equal to the given value)
    * right-tailed test (null hypothesis  -  the  median  is  less than or
      equal to the given value)

Requirements:
    * the scale of measurement should be ordinal, interval or  ratio (i.e.
      the test could not be applied to nominal variables).
    * the distribution should be continuous and symmetric relative to  its
      median.
    * number of distinct values in the X array should be greater than 4

The test is non-parametric and doesn't require distribution X to be normal

Input parameters:
    X       -   sample. Array whose index goes from 0 to N-1.
    N       -   size of the sample.
    Median  -   assumed median value.

Output parameters:
    BothTails   -   p-value for two-tailed test.
                    If BothTails is less than the given significance level
                    the null hypothesis is rejected.
    LeftTail    -   p-value for left-tailed test.
                    If LeftTail is less than the given significance level,
                    the null hypothesis is rejected.
    RightTail   -   p-value for right-tailed test.
                    If RightTail is less than the given significance level
                    the null hypothesis is rejected.

To calculate p-values, special approximation is used. This method lets  us
calculate p-values with two decimal places in interval [0.0001, 1].

"Two decimal places" does not sound very impressive, but in  practice  the
relative error of less than 1% is enough to make a decision.

There is no approximation outside the [0.0001, 1] interval. Therefore,  if
the significance level outlies this interval, the test returns 0.0001.

  -- ALGLIB --
     Copyright 08.09.2006 by Bochkanov Sergey
*************************************************************************/
void wilcoxonsignedranktest(const real_1d_array &x, const ae_int_t n, const double e, double &bothtails, double &lefttail, double &righttail);
}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (FUNCTIONS)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
void samplemoments(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* mean,
     double* variance,
     double* skewness,
     double* kurtosis,
     ae_state *_state);
double samplemean(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);
double samplevariance(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);
double sampleskewness(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);
double samplekurtosis(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);
void sampleadev(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* adev,
     ae_state *_state);
void samplemedian(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* median,
     ae_state *_state);
void samplepercentile(/* Real    */ ae_vector* x,
     ae_int_t n,
     double p,
     double* v,
     ae_state *_state);
double cov2(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_state *_state);
double pearsoncorr2(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_state *_state);
double spearmancorr2(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_state *_state);
void covm(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_matrix* c,
     ae_state *_state);
void _pexec_covm(/* Real    */ ae_matrix* x,
    ae_int_t n,
    ae_int_t m,
    /* Real    */ ae_matrix* c, ae_state *_state);
void pearsoncorrm(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_matrix* c,
     ae_state *_state);
void _pexec_pearsoncorrm(/* Real    */ ae_matrix* x,
    ae_int_t n,
    ae_int_t m,
    /* Real    */ ae_matrix* c, ae_state *_state);
void spearmancorrm(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_matrix* c,
     ae_state *_state);
void _pexec_spearmancorrm(/* Real    */ ae_matrix* x,
    ae_int_t n,
    ae_int_t m,
    /* Real    */ ae_matrix* c, ae_state *_state);
void covm2(/* Real    */ ae_matrix* x,
     /* Real    */ ae_matrix* y,
     ae_int_t n,
     ae_int_t m1,
     ae_int_t m2,
     /* Real    */ ae_matrix* c,
     ae_state *_state);
void _pexec_covm2(/* Real    */ ae_matrix* x,
    /* Real    */ ae_matrix* y,
    ae_int_t n,
    ae_int_t m1,
    ae_int_t m2,
    /* Real    */ ae_matrix* c, ae_state *_state);
void pearsoncorrm2(/* Real    */ ae_matrix* x,
     /* Real    */ ae_matrix* y,
     ae_int_t n,
     ae_int_t m1,
     ae_int_t m2,
     /* Real    */ ae_matrix* c,
     ae_state *_state);
void _pexec_pearsoncorrm2(/* Real    */ ae_matrix* x,
    /* Real    */ ae_matrix* y,
    ae_int_t n,
    ae_int_t m1,
    ae_int_t m2,
    /* Real    */ ae_matrix* c, ae_state *_state);
void spearmancorrm2(/* Real    */ ae_matrix* x,
     /* Real    */ ae_matrix* y,
     ae_int_t n,
     ae_int_t m1,
     ae_int_t m2,
     /* Real    */ ae_matrix* c,
     ae_state *_state);
void _pexec_spearmancorrm2(/* Real    */ ae_matrix* x,
    /* Real    */ ae_matrix* y,
    ae_int_t n,
    ae_int_t m1,
    ae_int_t m2,
    /* Real    */ ae_matrix* c, ae_state *_state);
void rankdata(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nfeatures,
     ae_state *_state);
void _pexec_rankdata(/* Real    */ ae_matrix* xy,
    ae_int_t npoints,
    ae_int_t nfeatures, ae_state *_state);
void rankdatacentered(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nfeatures,
     ae_state *_state);
void _pexec_rankdatacentered(/* Real    */ ae_matrix* xy,
    ae_int_t npoints,
    ae_int_t nfeatures, ae_state *_state);
double pearsoncorrelation(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_state *_state);
double spearmanrankcorrelation(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_state *_state);
void pearsoncorrelationsignificance(double r,
     ae_int_t n,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state);
void spearmanrankcorrelationsignificance(double r,
     ae_int_t n,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state);
void jarqueberatest(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* p,
     ae_state *_state);
void mannwhitneyutest(/* Real    */ ae_vector* x,
     ae_int_t n,
     /* Real    */ ae_vector* y,
     ae_int_t m,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state);
void onesamplesigntest(/* Real    */ ae_vector* x,
     ae_int_t n,
     double median,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state);
void studentttest1(/* Real    */ ae_vector* x,
     ae_int_t n,
     double mean,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state);
void studentttest2(/* Real    */ ae_vector* x,
     ae_int_t n,
     /* Real    */ ae_vector* y,
     ae_int_t m,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state);
void unequalvariancettest(/* Real    */ ae_vector* x,
     ae_int_t n,
     /* Real    */ ae_vector* y,
     ae_int_t m,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state);
void ftest(/* Real    */ ae_vector* x,
     ae_int_t n,
     /* Real    */ ae_vector* y,
     ae_int_t m,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state);
void onesamplevariancetest(/* Real    */ ae_vector* x,
     ae_int_t n,
     double variance,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state);
void wilcoxonsignedranktest(/* Real    */ ae_vector* x,
     ae_int_t n,
     double e,
     double* bothtails,
     double* lefttail,
     double* righttail,
     ae_state *_state);

}
#endif

