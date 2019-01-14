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
#ifndef _dataanalysis_pkg_h
#define _dataanalysis_pkg_h
#include "ap.h"
#include "alglibinternal.h"
#include "linalg.h"
#include "statistics.h"
#include "alglibmisc.h"
#include "specialfunctions.h"
#include "solvers.h"
#include "optimization.h"

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (DATATYPES)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
typedef struct
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
} cvreport;
typedef struct
{
    ae_matrix ct;
    ae_matrix ctbest;
    ae_vector xycbest;
    ae_vector xycprev;
    ae_vector d2;
    ae_vector csizes;
    apbuffers initbuf;
    ae_shared_pool updatepool;
} kmeansbuffers;
typedef struct
{
    ae_int_t npoints;
    ae_int_t nfeatures;
    ae_int_t disttype;
    ae_matrix xy;
    ae_matrix d;
    ae_int_t ahcalgo;
    ae_int_t kmeansrestarts;
    ae_int_t kmeansmaxits;
    ae_int_t kmeansinitalgo;
    ae_bool kmeansdbgnoits;
    ae_matrix tmpd;
    apbuffers distbuf;
    kmeansbuffers kmeanstmp;
} clusterizerstate;
typedef struct
{
    ae_int_t terminationtype;
    ae_int_t npoints;
    ae_vector p;
    ae_matrix z;
    ae_matrix pz;
    ae_matrix pm;
    ae_vector mergedist;
} ahcreport;
typedef struct
{
    ae_int_t npoints;
    ae_int_t nfeatures;
    ae_int_t terminationtype;
    ae_int_t iterationscount;
    double energy;
    ae_int_t k;
    ae_matrix c;
    ae_vector cidx;
} kmeansreport;
typedef struct
{
    ae_int_t nvars;
    ae_int_t nclasses;
    ae_int_t ntrees;
    ae_int_t bufsize;
    ae_vector trees;
} decisionforest;
typedef struct
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double oobrelclserror;
    double oobavgce;
    double oobrmserror;
    double oobavgerror;
    double oobavgrelerror;
} dfreport;
typedef struct
{
    ae_vector treebuf;
    ae_vector idxbuf;
    ae_vector tmpbufr;
    ae_vector tmpbufr2;
    ae_vector tmpbufi;
    ae_vector classibuf;
    ae_vector sortrbuf;
    ae_vector sortrbuf2;
    ae_vector sortibuf;
    ae_vector varpool;
    ae_vector evsbin;
    ae_vector evssplits;
} dfinternalbuffers;
typedef struct
{
    ae_vector w;
} linearmodel;
typedef struct
{
    ae_matrix c;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double cvrmserror;
    double cvavgerror;
    double cvavgrelerror;
    ae_int_t ncvdefects;
    ae_vector cvdefects;
} lrreport;
typedef struct
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
} modelerrors;
typedef struct
{
    double f;
    ae_vector g;
} smlpgrad;
typedef struct
{
    ae_int_t hlnetworktype;
    ae_int_t hlnormtype;
    ae_vector hllayersizes;
    ae_vector hlconnections;
    ae_vector hlneurons;
    ae_vector structinfo;
    ae_vector weights;
    ae_vector columnmeans;
    ae_vector columnsigmas;
    ae_vector neurons;
    ae_vector dfdnet;
    ae_vector derror;
    ae_vector x;
    ae_vector y;
    ae_matrix xy;
    ae_vector xyrow;
    ae_vector nwbuf;
    ae_vector integerbuf;
    modelerrors err;
    ae_vector rndbuf;
    ae_shared_pool buf;
    ae_shared_pool gradbuf;
    ae_matrix dummydxy;
    sparsematrix dummysxy;
    ae_vector dummyidx;
    ae_shared_pool dummypool;
} multilayerperceptron;
typedef struct
{
    ae_vector w;
} logitmodel;
typedef struct
{
    ae_bool brackt;
    ae_bool stage1;
    ae_int_t infoc;
    double dg;
    double dgm;
    double dginit;
    double dgtest;
    double dgx;
    double dgxm;
    double dgy;
    double dgym;
    double finit;
    double ftest1;
    double fm;
    double fx;
    double fxm;
    double fy;
    double fym;
    double stx;
    double sty;
    double stmin;
    double stmax;
    double width;
    double width1;
    double xtrapf;
} logitmcstate;
typedef struct
{
    ae_int_t ngrad;
    ae_int_t nhess;
} mnlreport;
typedef struct
{
    ae_int_t n;
    ae_vector states;
    ae_int_t npairs;
    ae_matrix data;
    ae_matrix ec;
    ae_matrix bndl;
    ae_matrix bndu;
    ae_matrix c;
    ae_vector ct;
    ae_int_t ccnt;
    ae_vector pw;
    ae_matrix priorp;
    double regterm;
    minbleicstate bs;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
    ae_int_t repnfev;
    ae_int_t repterminationtype;
    minbleicreport br;
    ae_vector tmpp;
    ae_vector effectivew;
    ae_vector effectivebndl;
    ae_vector effectivebndu;
    ae_matrix effectivec;
    ae_vector effectivect;
    ae_vector h;
    ae_matrix p;
} mcpdstate;
typedef struct
{
    ae_int_t inneriterationscount;
    ae_int_t outeriterationscount;
    ae_int_t nfev;
    ae_int_t terminationtype;
} mcpdreport;
typedef struct
{
    ae_int_t ensemblesize;
    ae_vector weights;
    ae_vector columnmeans;
    ae_vector columnsigmas;
    multilayerperceptron network;
    ae_vector y;
} mlpensemble;
typedef struct
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
    ae_int_t ngrad;
    ae_int_t nhess;
    ae_int_t ncholesky;
} mlpreport;
typedef struct
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
} mlpcvreport;
typedef struct
{
    ae_vector bestparameters;
    double bestrmserror;
    ae_bool randomizenetwork;
    multilayerperceptron network;
    minlbfgsstate optimizer;
    minlbfgsreport optimizerrep;
    ae_vector wbuf0;
    ae_vector wbuf1;
    ae_vector allminibatches;
    ae_vector currentminibatch;
    rcommstate rstate;
    ae_int_t algoused;
    ae_int_t minibatchsize;
    hqrndstate generator;
} smlptrnsession;
typedef struct
{
    ae_vector trnsubset;
    ae_vector valsubset;
    ae_shared_pool mlpsessions;
    mlpreport mlprep;
    multilayerperceptron network;
} mlpetrnsession;
typedef struct
{
    ae_int_t nin;
    ae_int_t nout;
    ae_bool rcpar;
    ae_int_t lbfgsfactor;
    double decay;
    double wstep;
    ae_int_t maxits;
    ae_int_t datatype;
    ae_int_t npoints;
    ae_matrix densexy;
    sparsematrix sparsexy;
    smlptrnsession session;
    ae_int_t ngradbatch;
    ae_vector subset;
    ae_int_t subsetsize;
    ae_vector valsubset;
    ae_int_t valsubsetsize;
    ae_int_t algokind;
    ae_int_t minibatchsize;
} mlptrainer;
typedef struct
{
    multilayerperceptron network;
    mlpreport rep;
    ae_vector subset;
    ae_int_t subsetsize;
    ae_vector xyrow;
    ae_vector y;
    ae_int_t ngrad;
    ae_shared_pool trnpool;
} mlpparallelizationcv;

}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS C++ INTERFACE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib
{



/*************************************************************************
This structure is a clusterization engine.

You should not try to access its fields directly.
Use ALGLIB functions in order to work with this object.

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
class _clusterizerstate_owner
{
public:
    _clusterizerstate_owner();
    _clusterizerstate_owner(const _clusterizerstate_owner &rhs);
    _clusterizerstate_owner& operator=(const _clusterizerstate_owner &rhs);
    virtual ~_clusterizerstate_owner();
    alglib_impl::clusterizerstate* c_ptr();
    alglib_impl::clusterizerstate* c_ptr() const;
protected:
    alglib_impl::clusterizerstate *p_struct;
};
class clusterizerstate : public _clusterizerstate_owner
{
public:
    clusterizerstate();
    clusterizerstate(const clusterizerstate &rhs);
    clusterizerstate& operator=(const clusterizerstate &rhs);
    virtual ~clusterizerstate();

};


/*************************************************************************
This structure  is used to store results of the agglomerative hierarchical
clustering (AHC).

Following information is returned:

* TerminationType - completion code:
  * 1   for successful completion of algorithm
  * -5  inappropriate combination of  clustering  algorithm  and  distance
        function was used. As for now, it  is  possible  only when  Ward's
        method is called for dataset with non-Euclidean distance function.
  In case negative completion code is returned,  other  fields  of  report
  structure are invalid and should not be used.

* NPoints contains number of points in the original dataset

* Z contains information about merges performed  (see below).  Z  contains
  indexes from the original (unsorted) dataset and it can be used when you
  need to know what points were merged. However, it is not convenient when
  you want to build a dendrograd (see below).

* if  you  want  to  build  dendrogram, you  can use Z, but it is not good
  option, because Z contains  indexes from  unsorted  dataset.  Dendrogram
  built from such dataset is likely to have intersections. So, you have to
  reorder you points before building dendrogram.
  Permutation which reorders point is returned in P. Another representation
  of  merges,  which  is  more  convenient for dendorgram construction, is
  returned in PM.

* more information on format of Z, P and PM can be found below and in the
  examples from ALGLIB Reference Manual.

FORMAL DESCRIPTION OF FIELDS:
    NPoints         number of points
    Z               array[NPoints-1,2],  contains   indexes   of  clusters
                    linked in pairs to  form  clustering  tree.  I-th  row
                    corresponds to I-th merge:
                    * Z[I,0] - index of the first cluster to merge
                    * Z[I,1] - index of the second cluster to merge
                    * Z[I,0]<Z[I,1]
                    * clusters are  numbered  from 0 to 2*NPoints-2,  with
                      indexes from 0 to NPoints-1 corresponding to  points
                      of the original dataset, and indexes from NPoints to
                      2*NPoints-2  correspond  to  clusters  generated  by
                      subsequent  merges  (I-th  row  of Z creates cluster
                      with index NPoints+I).

                    IMPORTANT: indexes in Z[] are indexes in the ORIGINAL,
                    unsorted dataset. In addition to  Z algorithm  outputs
                    permutation which rearranges points in such  way  that
                    subsequent merges are  performed  on  adjacent  points
                    (such order is needed if you want to build dendrogram).
                    However,  indexes  in  Z  are  related  to   original,
                    unrearranged sequence of points.

    P               array[NPoints], permutation which reorders points  for
                    dendrogram  construction.  P[i] contains  index of the
                    position  where  we  should  move  I-th  point  of the
                    original dataset in order to apply merges PZ/PM.

    PZ              same as Z, but for permutation of points given  by  P.
                    The  only  thing  which  changed  are  indexes  of the
                    original points; indexes of clusters remained same.

    MergeDist       array[NPoints-1], contains distances between  clusters
                    being merged (MergeDist[i] correspond to merge  stored
                    in Z[i,...]):
                    * CLINK, SLINK and  average  linkage algorithms report
                      "raw", unmodified distance metric.
                    * Ward's   method   reports   weighted   intra-cluster
                      variance, which is equal to ||Ca-Cb||^2 * Sa*Sb/(Sa+Sb).
                      Here  A  and  B  are  clusters being merged, Ca is a
                      center of A, Cb is a center of B, Sa is a size of A,
                      Sb is a size of B.

    PM              array[NPoints-1,6], another representation of  merges,
                    which is suited for dendrogram construction. It  deals
                    with rearranged points (permutation P is applied)  and
                    represents merges in a form which different  from  one
                    used by Z.
                    For each I from 0 to NPoints-2, I-th row of PM represents
                    merge performed on two clusters C0 and C1. Here:
                    * C0 contains points with indexes PM[I,0]...PM[I,1]
                    * C1 contains points with indexes PM[I,2]...PM[I,3]
                    * indexes stored in PM are given for dataset sorted
                      according to permutation P
                    * PM[I,1]=PM[I,2]-1 (only adjacent clusters are merged)
                    * PM[I,0]<=PM[I,1], PM[I,2]<=PM[I,3], i.e. both
                      clusters contain at least one point
                    * heights of "subdendrograms" corresponding  to  C0/C1
                      are stored in PM[I,4]  and  PM[I,5].  Subdendrograms
                      corresponding   to   single-point   clusters    have
                      height=0. Dendrogram of the merge result has  height
                      H=max(H0,H1)+1.

NOTE: there is one-to-one correspondence between merges described by Z and
      PM. I-th row of Z describes same merge of clusters as I-th row of PM,
      with "left" cluster from Z corresponding to the "left" one from PM.

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
class _ahcreport_owner
{
public:
    _ahcreport_owner();
    _ahcreport_owner(const _ahcreport_owner &rhs);
    _ahcreport_owner& operator=(const _ahcreport_owner &rhs);
    virtual ~_ahcreport_owner();
    alglib_impl::ahcreport* c_ptr();
    alglib_impl::ahcreport* c_ptr() const;
protected:
    alglib_impl::ahcreport *p_struct;
};
class ahcreport : public _ahcreport_owner
{
public:
    ahcreport();
    ahcreport(const ahcreport &rhs);
    ahcreport& operator=(const ahcreport &rhs);
    virtual ~ahcreport();
    ae_int_t &terminationtype;
    ae_int_t &npoints;
    integer_1d_array p;
    integer_2d_array z;
    integer_2d_array pz;
    integer_2d_array pm;
    real_1d_array mergedist;

};


/*************************************************************************
This  structure   is  used  to  store  results of the  k-means  clustering
algorithm.

Following information is always returned:
* NPoints contains number of points in the original dataset
* TerminationType contains completion code, negative on failure, positive
  on success
* K contains number of clusters

For positive TerminationType we return:
* NFeatures contains number of variables in the original dataset
* C, which contains centers found by algorithm
* CIdx, which maps points of the original dataset to clusters

FORMAL DESCRIPTION OF FIELDS:
    NPoints         number of points, >=0
    NFeatures       number of variables, >=1
    TerminationType completion code:
                    * -5 if  distance  type  is  anything  different  from
                         Euclidean metric
                    * -3 for degenerate dataset: a) less  than  K  distinct
                         points, b) K=0 for non-empty dataset.
                    * +1 for successful completion
    K               number of clusters
    C               array[K,NFeatures], rows of the array store centers
    CIdx            array[NPoints], which contains cluster indexes
    IterationsCount actual number of iterations performed by clusterizer.
                    If algorithm performed more than one random restart,
                    total number of iterations is returned.
    Energy          merit function, "energy", sum  of  squared  deviations
                    from cluster centers

  -- ALGLIB --
     Copyright 27.11.2012 by Bochkanov Sergey
*************************************************************************/
class _kmeansreport_owner
{
public:
    _kmeansreport_owner();
    _kmeansreport_owner(const _kmeansreport_owner &rhs);
    _kmeansreport_owner& operator=(const _kmeansreport_owner &rhs);
    virtual ~_kmeansreport_owner();
    alglib_impl::kmeansreport* c_ptr();
    alglib_impl::kmeansreport* c_ptr() const;
protected:
    alglib_impl::kmeansreport *p_struct;
};
class kmeansreport : public _kmeansreport_owner
{
public:
    kmeansreport();
    kmeansreport(const kmeansreport &rhs);
    kmeansreport& operator=(const kmeansreport &rhs);
    virtual ~kmeansreport();
    ae_int_t &npoints;
    ae_int_t &nfeatures;
    ae_int_t &terminationtype;
    ae_int_t &iterationscount;
    double &energy;
    ae_int_t &k;
    real_2d_array c;
    integer_1d_array cidx;

};



/*************************************************************************

*************************************************************************/
class _decisionforest_owner
{
public:
    _decisionforest_owner();
    _decisionforest_owner(const _decisionforest_owner &rhs);
    _decisionforest_owner& operator=(const _decisionforest_owner &rhs);
    virtual ~_decisionforest_owner();
    alglib_impl::decisionforest* c_ptr();
    alglib_impl::decisionforest* c_ptr() const;
protected:
    alglib_impl::decisionforest *p_struct;
};
class decisionforest : public _decisionforest_owner
{
public:
    decisionforest();
    decisionforest(const decisionforest &rhs);
    decisionforest& operator=(const decisionforest &rhs);
    virtual ~decisionforest();

};


/*************************************************************************

*************************************************************************/
class _dfreport_owner
{
public:
    _dfreport_owner();
    _dfreport_owner(const _dfreport_owner &rhs);
    _dfreport_owner& operator=(const _dfreport_owner &rhs);
    virtual ~_dfreport_owner();
    alglib_impl::dfreport* c_ptr();
    alglib_impl::dfreport* c_ptr() const;
protected:
    alglib_impl::dfreport *p_struct;
};
class dfreport : public _dfreport_owner
{
public:
    dfreport();
    dfreport(const dfreport &rhs);
    dfreport& operator=(const dfreport &rhs);
    virtual ~dfreport();
    double &relclserror;
    double &avgce;
    double &rmserror;
    double &avgerror;
    double &avgrelerror;
    double &oobrelclserror;
    double &oobavgce;
    double &oobrmserror;
    double &oobavgerror;
    double &oobavgrelerror;

};

/*************************************************************************

*************************************************************************/
class _linearmodel_owner
{
public:
    _linearmodel_owner();
    _linearmodel_owner(const _linearmodel_owner &rhs);
    _linearmodel_owner& operator=(const _linearmodel_owner &rhs);
    virtual ~_linearmodel_owner();
    alglib_impl::linearmodel* c_ptr();
    alglib_impl::linearmodel* c_ptr() const;
protected:
    alglib_impl::linearmodel *p_struct;
};
class linearmodel : public _linearmodel_owner
{
public:
    linearmodel();
    linearmodel(const linearmodel &rhs);
    linearmodel& operator=(const linearmodel &rhs);
    virtual ~linearmodel();

};


/*************************************************************************
LRReport structure contains additional information about linear model:
* C             -   covariation matrix,  array[0..NVars,0..NVars].
                    C[i,j] = Cov(A[i],A[j])
* RMSError      -   root mean square error on a training set
* AvgError      -   average error on a training set
* AvgRelError   -   average relative error on a training set (excluding
                    observations with zero function value).
* CVRMSError    -   leave-one-out cross-validation estimate of
                    generalization error. Calculated using fast algorithm
                    with O(NVars*NPoints) complexity.
* CVAvgError    -   cross-validation estimate of average error
* CVAvgRelError -   cross-validation estimate of average relative error

All other fields of the structure are intended for internal use and should
not be used outside ALGLIB.
*************************************************************************/
class _lrreport_owner
{
public:
    _lrreport_owner();
    _lrreport_owner(const _lrreport_owner &rhs);
    _lrreport_owner& operator=(const _lrreport_owner &rhs);
    virtual ~_lrreport_owner();
    alglib_impl::lrreport* c_ptr();
    alglib_impl::lrreport* c_ptr() const;
protected:
    alglib_impl::lrreport *p_struct;
};
class lrreport : public _lrreport_owner
{
public:
    lrreport();
    lrreport(const lrreport &rhs);
    lrreport& operator=(const lrreport &rhs);
    virtual ~lrreport();
    real_2d_array c;
    double &rmserror;
    double &avgerror;
    double &avgrelerror;
    double &cvrmserror;
    double &cvavgerror;
    double &cvavgrelerror;
    ae_int_t &ncvdefects;
    integer_1d_array cvdefects;

};





/*************************************************************************
Model's errors:
    * RelCLSError   -   fraction of misclassified cases.
    * AvgCE         -   acerage cross-entropy
    * RMSError      -   root-mean-square error
    * AvgError      -   average error
    * AvgRelError   -   average relative error

NOTE 1: RelCLSError/AvgCE are zero on regression problems.

NOTE 2: on classification problems  RMSError/AvgError/AvgRelError  contain
        errors in prediction of posterior probabilities
*************************************************************************/
class _modelerrors_owner
{
public:
    _modelerrors_owner();
    _modelerrors_owner(const _modelerrors_owner &rhs);
    _modelerrors_owner& operator=(const _modelerrors_owner &rhs);
    virtual ~_modelerrors_owner();
    alglib_impl::modelerrors* c_ptr();
    alglib_impl::modelerrors* c_ptr() const;
protected:
    alglib_impl::modelerrors *p_struct;
};
class modelerrors : public _modelerrors_owner
{
public:
    modelerrors();
    modelerrors(const modelerrors &rhs);
    modelerrors& operator=(const modelerrors &rhs);
    virtual ~modelerrors();
    double &relclserror;
    double &avgce;
    double &rmserror;
    double &avgerror;
    double &avgrelerror;

};


/*************************************************************************

*************************************************************************/
class _multilayerperceptron_owner
{
public:
    _multilayerperceptron_owner();
    _multilayerperceptron_owner(const _multilayerperceptron_owner &rhs);
    _multilayerperceptron_owner& operator=(const _multilayerperceptron_owner &rhs);
    virtual ~_multilayerperceptron_owner();
    alglib_impl::multilayerperceptron* c_ptr();
    alglib_impl::multilayerperceptron* c_ptr() const;
protected:
    alglib_impl::multilayerperceptron *p_struct;
};
class multilayerperceptron : public _multilayerperceptron_owner
{
public:
    multilayerperceptron();
    multilayerperceptron(const multilayerperceptron &rhs);
    multilayerperceptron& operator=(const multilayerperceptron &rhs);
    virtual ~multilayerperceptron();

};

/*************************************************************************

*************************************************************************/
class _logitmodel_owner
{
public:
    _logitmodel_owner();
    _logitmodel_owner(const _logitmodel_owner &rhs);
    _logitmodel_owner& operator=(const _logitmodel_owner &rhs);
    virtual ~_logitmodel_owner();
    alglib_impl::logitmodel* c_ptr();
    alglib_impl::logitmodel* c_ptr() const;
protected:
    alglib_impl::logitmodel *p_struct;
};
class logitmodel : public _logitmodel_owner
{
public:
    logitmodel();
    logitmodel(const logitmodel &rhs);
    logitmodel& operator=(const logitmodel &rhs);
    virtual ~logitmodel();

};


/*************************************************************************
MNLReport structure contains information about training process:
* NGrad     -   number of gradient calculations
* NHess     -   number of Hessian calculations
*************************************************************************/
class _mnlreport_owner
{
public:
    _mnlreport_owner();
    _mnlreport_owner(const _mnlreport_owner &rhs);
    _mnlreport_owner& operator=(const _mnlreport_owner &rhs);
    virtual ~_mnlreport_owner();
    alglib_impl::mnlreport* c_ptr();
    alglib_impl::mnlreport* c_ptr() const;
protected:
    alglib_impl::mnlreport *p_struct;
};
class mnlreport : public _mnlreport_owner
{
public:
    mnlreport();
    mnlreport(const mnlreport &rhs);
    mnlreport& operator=(const mnlreport &rhs);
    virtual ~mnlreport();
    ae_int_t &ngrad;
    ae_int_t &nhess;

};

/*************************************************************************
This structure is a MCPD (Markov Chains for Population Data) solver.

You should use ALGLIB functions in order to work with this object.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
class _mcpdstate_owner
{
public:
    _mcpdstate_owner();
    _mcpdstate_owner(const _mcpdstate_owner &rhs);
    _mcpdstate_owner& operator=(const _mcpdstate_owner &rhs);
    virtual ~_mcpdstate_owner();
    alglib_impl::mcpdstate* c_ptr();
    alglib_impl::mcpdstate* c_ptr() const;
protected:
    alglib_impl::mcpdstate *p_struct;
};
class mcpdstate : public _mcpdstate_owner
{
public:
    mcpdstate();
    mcpdstate(const mcpdstate &rhs);
    mcpdstate& operator=(const mcpdstate &rhs);
    virtual ~mcpdstate();

};


/*************************************************************************
This structure is a MCPD training report:
    InnerIterationsCount    -   number of inner iterations of the
                                underlying optimization algorithm
    OuterIterationsCount    -   number of outer iterations of the
                                underlying optimization algorithm
    NFEV                    -   number of merit function evaluations
    TerminationType         -   termination type
                                (same as for MinBLEIC optimizer, positive
                                values denote success, negative ones -
                                failure)

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
class _mcpdreport_owner
{
public:
    _mcpdreport_owner();
    _mcpdreport_owner(const _mcpdreport_owner &rhs);
    _mcpdreport_owner& operator=(const _mcpdreport_owner &rhs);
    virtual ~_mcpdreport_owner();
    alglib_impl::mcpdreport* c_ptr();
    alglib_impl::mcpdreport* c_ptr() const;
protected:
    alglib_impl::mcpdreport *p_struct;
};
class mcpdreport : public _mcpdreport_owner
{
public:
    mcpdreport();
    mcpdreport(const mcpdreport &rhs);
    mcpdreport& operator=(const mcpdreport &rhs);
    virtual ~mcpdreport();
    ae_int_t &inneriterationscount;
    ae_int_t &outeriterationscount;
    ae_int_t &nfev;
    ae_int_t &terminationtype;

};

/*************************************************************************
Neural networks ensemble
*************************************************************************/
class _mlpensemble_owner
{
public:
    _mlpensemble_owner();
    _mlpensemble_owner(const _mlpensemble_owner &rhs);
    _mlpensemble_owner& operator=(const _mlpensemble_owner &rhs);
    virtual ~_mlpensemble_owner();
    alglib_impl::mlpensemble* c_ptr();
    alglib_impl::mlpensemble* c_ptr() const;
protected:
    alglib_impl::mlpensemble *p_struct;
};
class mlpensemble : public _mlpensemble_owner
{
public:
    mlpensemble();
    mlpensemble(const mlpensemble &rhs);
    mlpensemble& operator=(const mlpensemble &rhs);
    virtual ~mlpensemble();

};

/*************************************************************************
Training report:
    * RelCLSError   -   fraction of misclassified cases.
    * AvgCE         -   acerage cross-entropy
    * RMSError      -   root-mean-square error
    * AvgError      -   average error
    * AvgRelError   -   average relative error
    * NGrad         -   number of gradient calculations
    * NHess         -   number of Hessian calculations
    * NCholesky     -   number of Cholesky decompositions

NOTE 1: RelCLSError/AvgCE are zero on regression problems.

NOTE 2: on classification problems  RMSError/AvgError/AvgRelError  contain
        errors in prediction of posterior probabilities
*************************************************************************/
class _mlpreport_owner
{
public:
    _mlpreport_owner();
    _mlpreport_owner(const _mlpreport_owner &rhs);
    _mlpreport_owner& operator=(const _mlpreport_owner &rhs);
    virtual ~_mlpreport_owner();
    alglib_impl::mlpreport* c_ptr();
    alglib_impl::mlpreport* c_ptr() const;
protected:
    alglib_impl::mlpreport *p_struct;
};
class mlpreport : public _mlpreport_owner
{
public:
    mlpreport();
    mlpreport(const mlpreport &rhs);
    mlpreport& operator=(const mlpreport &rhs);
    virtual ~mlpreport();
    double &relclserror;
    double &avgce;
    double &rmserror;
    double &avgerror;
    double &avgrelerror;
    ae_int_t &ngrad;
    ae_int_t &nhess;
    ae_int_t &ncholesky;

};


/*************************************************************************
Cross-validation estimates of generalization error
*************************************************************************/
class _mlpcvreport_owner
{
public:
    _mlpcvreport_owner();
    _mlpcvreport_owner(const _mlpcvreport_owner &rhs);
    _mlpcvreport_owner& operator=(const _mlpcvreport_owner &rhs);
    virtual ~_mlpcvreport_owner();
    alglib_impl::mlpcvreport* c_ptr();
    alglib_impl::mlpcvreport* c_ptr() const;
protected:
    alglib_impl::mlpcvreport *p_struct;
};
class mlpcvreport : public _mlpcvreport_owner
{
public:
    mlpcvreport();
    mlpcvreport(const mlpcvreport &rhs);
    mlpcvreport& operator=(const mlpcvreport &rhs);
    virtual ~mlpcvreport();
    double &relclserror;
    double &avgce;
    double &rmserror;
    double &avgerror;
    double &avgrelerror;

};


/*************************************************************************
Trainer object for neural network.

You should not try to access fields of this object directly -  use  ALGLIB
functions to work with this object.
*************************************************************************/
class _mlptrainer_owner
{
public:
    _mlptrainer_owner();
    _mlptrainer_owner(const _mlptrainer_owner &rhs);
    _mlptrainer_owner& operator=(const _mlptrainer_owner &rhs);
    virtual ~_mlptrainer_owner();
    alglib_impl::mlptrainer* c_ptr();
    alglib_impl::mlptrainer* c_ptr() const;
protected:
    alglib_impl::mlptrainer *p_struct;
};
class mlptrainer : public _mlptrainer_owner
{
public:
    mlptrainer();
    mlptrainer(const mlptrainer &rhs);
    mlptrainer& operator=(const mlptrainer &rhs);
    virtual ~mlptrainer();

};

/*************************************************************************
Optimal binary classification

Algorithms finds optimal (=with minimal cross-entropy) binary partition.
Internal subroutine.

INPUT PARAMETERS:
    A       -   array[0..N-1], variable
    C       -   array[0..N-1], class numbers (0 or 1).
    N       -   array size

OUTPUT PARAMETERS:
    Info    -   completetion code:
                * -3, all values of A[] are same (partition is impossible)
                * -2, one of C[] is incorrect (<0, >1)
                * -1, incorrect pararemets were passed (N<=0).
                *  1, OK
    Threshold-  partiton boundary. Left part contains values which are
                strictly less than Threshold. Right part contains values
                which are greater than or equal to Threshold.
    PAL, PBL-   probabilities P(0|v<Threshold) and P(1|v<Threshold)
    PAR, PBR-   probabilities P(0|v>=Threshold) and P(1|v>=Threshold)
    CVE     -   cross-validation estimate of cross-entropy

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsoptimalsplit2(const real_1d_array &a, const integer_1d_array &c, const ae_int_t n, ae_int_t &info, double &threshold, double &pal, double &pbl, double &par, double &pbr, double &cve);


/*************************************************************************
Optimal partition, internal subroutine. Fast version.

Accepts:
    A       array[0..N-1]       array of attributes     array[0..N-1]
    C       array[0..N-1]       array of class labels
    TiesBuf array[0..N]         temporaries (ties)
    CntBuf  array[0..2*NC-1]    temporaries (counts)
    Alpha                       centering factor (0<=alpha<=1, recommended value - 0.05)
    BufR    array[0..N-1]       temporaries
    BufI    array[0..N-1]       temporaries

Output:
    Info    error code (">0"=OK, "<0"=bad)
    RMS     training set RMS error
    CVRMS   leave-one-out RMS error

Note:
    content of all arrays is changed by subroutine;
    it doesn't allocate temporaries.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void dsoptimalsplit2fast(real_1d_array &a, integer_1d_array &c, integer_1d_array &tiesbuf, integer_1d_array &cntbuf, real_1d_array &bufr, integer_1d_array &bufi, const ae_int_t n, const ae_int_t nc, const double alpha, ae_int_t &info, double &threshold, double &rms, double &cvrms);

/*************************************************************************
This function initializes clusterizer object. Newly initialized object  is
empty, i.e. it does not contain dataset. You should use it as follows:
1. creation
2. dataset is added with ClusterizerSetPoints()
3. additional parameters are set
3. clusterization is performed with one of the clustering functions

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizercreate(clusterizerstate &s);


/*************************************************************************
This function adds dataset to the clusterizer structure.

This function overrides all previous calls  of  ClusterizerSetPoints()  or
ClusterizerSetDistances().

INPUT PARAMETERS:
    S       -   clusterizer state, initialized by ClusterizerCreate()
    XY      -   array[NPoints,NFeatures], dataset
    NPoints -   number of points, >=0
    NFeatures-  number of features, >=1
    DistType-   distance function:
                *  0    Chebyshev distance  (L-inf norm)
                *  1    city block distance (L1 norm)
                *  2    Euclidean distance  (L2 norm), non-squared
                * 10    Pearson correlation:
                        dist(a,b) = 1-corr(a,b)
                * 11    Absolute Pearson correlation:
                        dist(a,b) = 1-|corr(a,b)|
                * 12    Uncentered Pearson correlation (cosine of the angle):
                        dist(a,b) = a'*b/(|a|*|b|)
                * 13    Absolute uncentered Pearson correlation
                        dist(a,b) = |a'*b|/(|a|*|b|)
                * 20    Spearman rank correlation:
                        dist(a,b) = 1-rankcorr(a,b)
                * 21    Absolute Spearman rank correlation
                        dist(a,b) = 1-|rankcorr(a,b)|

NOTE 1: different distance functions have different performance penalty:
        * Euclidean or Pearson correlation distances are the fastest ones
        * Spearman correlation distance function is a bit slower
        * city block and Chebyshev distances are order of magnitude slower

        The reason behing difference in performance is that correlation-based
        distance functions are computed using optimized linear algebra kernels,
        while Chebyshev and city block distance functions are computed using
        simple nested loops with two branches at each iteration.

NOTE 2: different clustering algorithms have different limitations:
        * agglomerative hierarchical clustering algorithms may be used with
          any kind of distance metric
        * k-means++ clustering algorithm may be used only  with  Euclidean
          distance function
        Thus, list of specific clustering algorithms you may  use  depends
        on distance function you specify when you set your dataset.

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizersetpoints(const clusterizerstate &s, const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype);
void clusterizersetpoints(const clusterizerstate &s, const real_2d_array &xy, const ae_int_t disttype);


/*************************************************************************
This function adds dataset given by distance  matrix  to  the  clusterizer
structure. It is important that dataset is not  given  explicitly  -  only
distance matrix is given.

This function overrides all previous calls  of  ClusterizerSetPoints()  or
ClusterizerSetDistances().

INPUT PARAMETERS:
    S       -   clusterizer state, initialized by ClusterizerCreate()
    D       -   array[NPoints,NPoints], distance matrix given by its upper
                or lower triangle (main diagonal is  ignored  because  its
                entries are expected to be zero).
    NPoints -   number of points
    IsUpper -   whether upper or lower triangle of D is given.

NOTE 1: different clustering algorithms have different limitations:
        * agglomerative hierarchical clustering algorithms may be used with
          any kind of distance metric, including one  which  is  given  by
          distance matrix
        * k-means++ clustering algorithm may be used only  with  Euclidean
          distance function and explicitly given points - it  can  not  be
          used with dataset given by distance matrix
        Thus, if you call this function, you will be unable to use k-means
        clustering algorithm to process your problem.

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizersetdistances(const clusterizerstate &s, const real_2d_array &d, const ae_int_t npoints, const bool isupper);
void clusterizersetdistances(const clusterizerstate &s, const real_2d_array &d, const bool isupper);


/*************************************************************************
This function sets agglomerative hierarchical clustering algorithm

INPUT PARAMETERS:
    S       -   clusterizer state, initialized by ClusterizerCreate()
    Algo    -   algorithm type:
                * 0     complete linkage (default algorithm)
                * 1     single linkage
                * 2     unweighted average linkage
                * 3     weighted average linkage
                * 4     Ward's method

NOTE: Ward's method works correctly only with Euclidean  distance,  that's
      why algorithm will return negative termination  code  (failure)  for
      any other distance type.

      It is possible, however,  to  use  this  method  with  user-supplied
      distance matrix. It  is  your  responsibility  to pass one which was
      calculated with Euclidean distance function.

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizersetahcalgo(const clusterizerstate &s, const ae_int_t algo);


/*************************************************************************
This  function  sets k-means properties:  number  of  restarts and maximum
number of iterations per one run.

INPUT PARAMETERS:
    S       -   clusterizer state, initialized by ClusterizerCreate()
    Restarts-   restarts count, >=1.
                k-means++ algorithm performs several restarts and  chooses
                best set of centers (one with minimum squared distance).
    MaxIts  -   maximum number of k-means iterations performed during  one
                run. >=0, zero value means that algorithm performs unlimited
                number of iterations.

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizersetkmeanslimits(const clusterizerstate &s, const ae_int_t restarts, const ae_int_t maxits);


/*************************************************************************
This function sets k-means  initialization  algorithm.  Several  different
algorithms can be chosen, including k-means++.

INPUT PARAMETERS:
    S       -   clusterizer state, initialized by ClusterizerCreate()
    InitAlgo-   initialization algorithm:
                * 0  automatic selection ( different  versions  of  ALGLIB
                     may select different algorithms)
                * 1  random initialization
                * 2  k-means++ initialization  (best  quality  of  initial
                     centers, but long  non-parallelizable  initialization
                     phase with bad cache locality)
                * 3  "fast-greedy"  algorithm  with  efficient,  easy   to
                     parallelize initialization. Quality of initial centers
                     is  somewhat  worse  than  that  of  k-means++.  This
                     algorithm is a default one in the current version  of
                     ALGLIB.
                *-1  "debug" algorithm which always selects first  K  rows
                     of dataset; this algorithm is used for debug purposes
                     only. Do not use it in the industrial code!

  -- ALGLIB --
     Copyright 21.01.2015 by Bochkanov Sergey
*************************************************************************/
void clusterizersetkmeansinit(const clusterizerstate &s, const ae_int_t initalgo);


/*************************************************************************
This function performs agglomerative hierarchical clustering

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Agglomerative  hierarchical  clustering  algorithm  has  two   phases:
  ! distance matrix calculation  and  clustering  itself. Only first phase
  ! (distance matrix calculation) is accelerated by Intel MKL  and  multi-
  ! threading. Thus, acceleration is significant only for  medium or high-
  ! dimensional problems.
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS:
    S       -   clusterizer state, initialized by ClusterizerCreate()

OUTPUT PARAMETERS:
    Rep     -   clustering results; see description of AHCReport
                structure for more information.

NOTE 1: hierarchical clustering algorithms require large amounts of memory.
        In particular, this implementation needs  sizeof(double)*NPoints^2
        bytes, which are used to store distance matrix. In  case  we  work
        with user-supplied matrix, this amount is multiplied by 2 (we have
        to store original matrix and to work with its copy).

        For example, problem with 10000 points  would require 800M of RAM,
        even when working in a 1-dimensional space.

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizerrunahc(const clusterizerstate &s, ahcreport &rep);
void smp_clusterizerrunahc(const clusterizerstate &s, ahcreport &rep);


/*************************************************************************
This function performs clustering by k-means++ algorithm.

You may change algorithm properties by calling:
* ClusterizerSetKMeansLimits() to change number of restarts or iterations
* ClusterizerSetKMeansInit() to change initialization algorithm

By  default,  one  restart  and  unlimited number of iterations are  used.
Initialization algorithm is chosen automatically.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes  two important  improvements  of
  ! this function:
  ! * multicore support (can be used from C# and C++)
  ! * access to high-performance C++ core (actual for C# users)
  !
  ! K-means clustering  algorithm has two  phases:  selection  of  initial
  ! centers  and  clustering  itself.  ALGLIB  parallelizes  both  phases.
  ! Parallel version is optimized for the following  scenario:  medium  or
  ! high-dimensional problem (20 or more dimensions) with large number  of
  ! points and clusters. However, some speed-up can be obtained even  when
  ! assumptions above are violated.
  !
  ! As for native-vs-managed comparison, working with native  core  brings
  ! 30-40% improvement in speed over pure C# version of ALGLIB.
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS:
    S       -   clusterizer state, initialized by ClusterizerCreate()
    K       -   number of clusters, K>=0.
                K  can  be  zero only when algorithm is called  for  empty
                dataset,  in   this   case   completion  code  is  set  to
                success (+1).
                If  K=0  and  dataset  size  is  non-zero,  we   can   not
                meaningfully assign points to some center  (there  are  no
                centers because K=0) and  return  -3  as  completion  code
                (failure).

OUTPUT PARAMETERS:
    Rep     -   clustering results; see description of KMeansReport
                structure for more information.

NOTE 1: k-means  clustering  can  be  performed  only  for  datasets  with
        Euclidean  distance  function.  Algorithm  will  return   negative
        completion code in Rep.TerminationType in case dataset  was  added
        to clusterizer with DistType other than Euclidean (or dataset  was
        specified by distance matrix instead of explicitly given points).

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizerrunkmeans(const clusterizerstate &s, const ae_int_t k, kmeansreport &rep);
void smp_clusterizerrunkmeans(const clusterizerstate &s, const ae_int_t k, kmeansreport &rep);


/*************************************************************************
This function returns distance matrix for dataset

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multicore support
  !
  ! Agglomerative  hierarchical  clustering  algorithm  has  two   phases:
  ! distance matrix calculation  and  clustering  itself. Only first phase
  ! (distance matrix calculation) is accelerated by Intel MKL  and  multi-
  ! threading. Thus, acceleration is significant only for  medium or high-
  ! dimensional problems.
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS:
    XY      -   array[NPoints,NFeatures], dataset
    NPoints -   number of points, >=0
    NFeatures-  number of features, >=1
    DistType-   distance function:
                *  0    Chebyshev distance  (L-inf norm)
                *  1    city block distance (L1 norm)
                *  2    Euclidean distance  (L2 norm, non-squared)
                * 10    Pearson correlation:
                        dist(a,b) = 1-corr(a,b)
                * 11    Absolute Pearson correlation:
                        dist(a,b) = 1-|corr(a,b)|
                * 12    Uncentered Pearson correlation (cosine of the angle):
                        dist(a,b) = a'*b/(|a|*|b|)
                * 13    Absolute uncentered Pearson correlation
                        dist(a,b) = |a'*b|/(|a|*|b|)
                * 20    Spearman rank correlation:
                        dist(a,b) = 1-rankcorr(a,b)
                * 21    Absolute Spearman rank correlation
                        dist(a,b) = 1-|rankcorr(a,b)|

OUTPUT PARAMETERS:
    D       -   array[NPoints,NPoints], distance matrix
                (full matrix is returned, with lower and upper triangles)

NOTE:  different distance functions have different performance penalty:
       * Euclidean or Pearson correlation distances are the fastest ones
       * Spearman correlation distance function is a bit slower
       * city block and Chebyshev distances are order of magnitude slower

       The reason behing difference in performance is that correlation-based
       distance functions are computed using optimized linear algebra kernels,
       while Chebyshev and city block distance functions are computed using
       simple nested loops with two branches at each iteration.

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizergetdistances(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype, real_2d_array &d);
void smp_clusterizergetdistances(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nfeatures, const ae_int_t disttype, real_2d_array &d);


/*************************************************************************
This function takes as input clusterization report Rep,  desired  clusters
count K, and builds top K clusters from hierarchical clusterization  tree.
It returns assignment of points to clusters (array of cluster indexes).

INPUT PARAMETERS:
    Rep     -   report from ClusterizerRunAHC() performed on XY
    K       -   desired number of clusters, 1<=K<=NPoints.
                K can be zero only when NPoints=0.

OUTPUT PARAMETERS:
    CIdx    -   array[NPoints], I-th element contains cluster index  (from
                0 to K-1) for I-th point of the dataset.
    CZ      -   array[K]. This array allows  to  convert  cluster  indexes
                returned by this function to indexes used by  Rep.Z.  J-th
                cluster returned by this function corresponds to  CZ[J]-th
                cluster stored in Rep.Z/PZ/PM.
                It is guaranteed that CZ[I]<CZ[I+1].

NOTE: K clusters built by this subroutine are assumed to have no hierarchy.
      Although  they  were  obtained  by  manipulation with top K nodes of
      dendrogram  (i.e.  hierarchical  decomposition  of  dataset),   this
      function does not return information about hierarchy.  Each  of  the
      clusters stand on its own.

NOTE: Cluster indexes returned by this function  does  not  correspond  to
      indexes returned in Rep.Z/PZ/PM. Either you work  with  hierarchical
      representation of the dataset (dendrogram), or you work with  "flat"
      representation returned by this function.  Each  of  representations
      has its own clusters indexing system (former uses [0, 2*NPoints-2]),
      while latter uses [0..K-1]), although  it  is  possible  to  perform
      conversion from one system to another by means of CZ array, returned
      by this function, which allows you to convert indexes stored in CIdx
      to the numeration system used by Rep.Z.

NOTE: this subroutine is optimized for moderate values of K. Say, for  K=5
      it will perform many times faster than  for  K=100.  Its  worst-case
      performance is O(N*K), although in average case  it  perform  better
      (up to O(N*log(K))).

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizergetkclusters(const ahcreport &rep, const ae_int_t k, integer_1d_array &cidx, integer_1d_array &cz);


/*************************************************************************
This  function  accepts  AHC  report  Rep,  desired  minimum  intercluster
distance and returns top clusters from  hierarchical  clusterization  tree
which are separated by distance R or HIGHER.

It returns assignment of points to clusters (array of cluster indexes).

There is one more function with similar name - ClusterizerSeparatedByCorr,
which returns clusters with intercluster correlation equal to R  or  LOWER
(note: higher for distance, lower for correlation).

INPUT PARAMETERS:
    Rep     -   report from ClusterizerRunAHC() performed on XY
    R       -   desired minimum intercluster distance, R>=0

OUTPUT PARAMETERS:
    K       -   number of clusters, 1<=K<=NPoints
    CIdx    -   array[NPoints], I-th element contains cluster index  (from
                0 to K-1) for I-th point of the dataset.
    CZ      -   array[K]. This array allows  to  convert  cluster  indexes
                returned by this function to indexes used by  Rep.Z.  J-th
                cluster returned by this function corresponds to  CZ[J]-th
                cluster stored in Rep.Z/PZ/PM.
                It is guaranteed that CZ[I]<CZ[I+1].

NOTE: K clusters built by this subroutine are assumed to have no hierarchy.
      Although  they  were  obtained  by  manipulation with top K nodes of
      dendrogram  (i.e.  hierarchical  decomposition  of  dataset),   this
      function does not return information about hierarchy.  Each  of  the
      clusters stand on its own.

NOTE: Cluster indexes returned by this function  does  not  correspond  to
      indexes returned in Rep.Z/PZ/PM. Either you work  with  hierarchical
      representation of the dataset (dendrogram), or you work with  "flat"
      representation returned by this function.  Each  of  representations
      has its own clusters indexing system (former uses [0, 2*NPoints-2]),
      while latter uses [0..K-1]), although  it  is  possible  to  perform
      conversion from one system to another by means of CZ array, returned
      by this function, which allows you to convert indexes stored in CIdx
      to the numeration system used by Rep.Z.

NOTE: this subroutine is optimized for moderate values of K. Say, for  K=5
      it will perform many times faster than  for  K=100.  Its  worst-case
      performance is O(N*K), although in average case  it  perform  better
      (up to O(N*log(K))).

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizerseparatedbydist(const ahcreport &rep, const double r, ae_int_t &k, integer_1d_array &cidx, integer_1d_array &cz);


/*************************************************************************
This  function  accepts  AHC  report  Rep,  desired  maximum  intercluster
correlation and returns top clusters from hierarchical clusterization tree
which are separated by correlation R or LOWER.

It returns assignment of points to clusters (array of cluster indexes).

There is one more function with similar name - ClusterizerSeparatedByDist,
which returns clusters with intercluster distance equal  to  R  or  HIGHER
(note: higher for distance, lower for correlation).

INPUT PARAMETERS:
    Rep     -   report from ClusterizerRunAHC() performed on XY
    R       -   desired maximum intercluster correlation, -1<=R<=+1

OUTPUT PARAMETERS:
    K       -   number of clusters, 1<=K<=NPoints
    CIdx    -   array[NPoints], I-th element contains cluster index  (from
                0 to K-1) for I-th point of the dataset.
    CZ      -   array[K]. This array allows  to  convert  cluster  indexes
                returned by this function to indexes used by  Rep.Z.  J-th
                cluster returned by this function corresponds to  CZ[J]-th
                cluster stored in Rep.Z/PZ/PM.
                It is guaranteed that CZ[I]<CZ[I+1].

NOTE: K clusters built by this subroutine are assumed to have no hierarchy.
      Although  they  were  obtained  by  manipulation with top K nodes of
      dendrogram  (i.e.  hierarchical  decomposition  of  dataset),   this
      function does not return information about hierarchy.  Each  of  the
      clusters stand on its own.

NOTE: Cluster indexes returned by this function  does  not  correspond  to
      indexes returned in Rep.Z/PZ/PM. Either you work  with  hierarchical
      representation of the dataset (dendrogram), or you work with  "flat"
      representation returned by this function.  Each  of  representations
      has its own clusters indexing system (former uses [0, 2*NPoints-2]),
      while latter uses [0..K-1]), although  it  is  possible  to  perform
      conversion from one system to another by means of CZ array, returned
      by this function, which allows you to convert indexes stored in CIdx
      to the numeration system used by Rep.Z.

NOTE: this subroutine is optimized for moderate values of K. Say, for  K=5
      it will perform many times faster than  for  K=100.  Its  worst-case
      performance is O(N*K), although in average case  it  perform  better
      (up to O(N*log(K))).

  -- ALGLIB --
     Copyright 10.07.2012 by Bochkanov Sergey
*************************************************************************/
void clusterizerseparatedbycorr(const ahcreport &rep, const double r, ae_int_t &k, integer_1d_array &cidx, integer_1d_array &cz);

/*************************************************************************
k-means++ clusterization.
Backward compatibility function, we recommend to use CLUSTERING subpackage
as better replacement.

  -- ALGLIB --
     Copyright 21.03.2009 by Bochkanov Sergey
*************************************************************************/
void kmeansgenerate(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t k, const ae_int_t restarts, ae_int_t &info, real_2d_array &c, integer_1d_array &xyc);

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
void dfserialize(decisionforest &obj, std::string &s_out);


/*************************************************************************
This function unserializes data structure from string.
*************************************************************************/
void dfunserialize(std::string &s_in, decisionforest &obj);


/*************************************************************************
This subroutine builds random decision forest.

INPUT PARAMETERS:
    XY          -   training set
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   task type:
                    * NClasses=1 - regression task with one
                                   dependent variable
                    * NClasses>1 - classification task with
                                   NClasses classes.
    NTrees      -   number of trees in a forest, NTrees>=1.
                    recommended values: 50-100.
    R           -   percent of a training set used to build
                    individual trees. 0<R<=1.
                    recommended values: 0.1 <= R <= 0.66.

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<1, NVars<1, NClasses<1, NTrees<1, R<=0
                          or R>1).
                    *  1, if task has been solved
    DF          -   model built
    Rep         -   training report, contains error on a training set
                    and out-of-bag estimates of generalization error.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfbuildrandomdecisionforest(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const double r, ae_int_t &info, decisionforest &df, dfreport &rep);


/*************************************************************************
This subroutine builds random decision forest.
This function gives ability to tune number of variables used when choosing
best split.

INPUT PARAMETERS:
    XY          -   training set
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   task type:
                    * NClasses=1 - regression task with one
                                   dependent variable
                    * NClasses>1 - classification task with
                                   NClasses classes.
    NTrees      -   number of trees in a forest, NTrees>=1.
                    recommended values: 50-100.
    NRndVars    -   number of variables used when choosing best split
    R           -   percent of a training set used to build
                    individual trees. 0<R<=1.
                    recommended values: 0.1 <= R <= 0.66.

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<1, NVars<1, NClasses<1, NTrees<1, R<=0
                          or R>1).
                    *  1, if task has been solved
    DF          -   model built
    Rep         -   training report, contains error on a training set
                    and out-of-bag estimates of generalization error.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfbuildrandomdecisionforestx1(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const ae_int_t nrndvars, const double r, ae_int_t &info, decisionforest &df, dfreport &rep);


/*************************************************************************
Procesing

INPUT PARAMETERS:
    DF      -   decision forest model
    X       -   input vector,  array[0..NVars-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

See also DFProcessI.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfprocess(const decisionforest &df, const real_1d_array &x, real_1d_array &y);


/*************************************************************************
'interactive' variant of DFProcess for languages like Python which support
constructs like "Y = DFProcessI(DF,X)" and interactive mode of interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void dfprocessi(const decisionforest &df, const real_1d_array &x, real_1d_array &y);


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfrelclserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgce(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for
    classification task, RMS error means error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfrmserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average relative error when estimating
    posterior probability of belonging to the correct class.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgrelerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);

/*************************************************************************
Linear regression

Subroutine builds model:

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1] + A(N)

and model found in ALGLIB format, covariation matrix, training set  errors
(rms,  average,  average  relative)   and  leave-one-out  cross-validation
estimate of the generalization error. CV  estimate calculated  using  fast
algorithm with O(NPoints*NVars) complexity.

When  covariation  matrix  is  calculated  standard deviations of function
values are assumed to be equal to RMS error on the training set.

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuild(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);


/*************************************************************************
Linear regression

Variant of LRBuild which uses vector of standatd deviations (errors in
function values).

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    S           -   standard deviations (errors in function values)
                    array[0..NPoints-1], S[i]>0.
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    * -2, if S[I]<=0
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuilds(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);


/*************************************************************************
Like LRBuildS, but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildzs(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);


/*************************************************************************
Like LRBuild but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildz(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);


/*************************************************************************
Unpacks coefficients of linear model.

INPUT PARAMETERS:
    LM          -   linear model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
                    constant term (intercept) is stored in the V[NVars].
    NVars       -   number of independent variables (one less than number
                    of coefficients)

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrunpack(const linearmodel &lm, real_1d_array &v, ae_int_t &nvars);


/*************************************************************************
"Packs" coefficients and creates linear model in ALGLIB format (LRUnpack
reversed).

INPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
    NVars       -   number of independent variables

OUTPUT PAREMETERS:
    LM          -   linear model.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrpack(const real_1d_array &v, const ae_int_t nvars, linearmodel &lm);


/*************************************************************************
Procesing

INPUT PARAMETERS:
    LM      -   linear model
    X       -   input vector,  array[0..NVars-1].

Result:
    value of linear model regression estimate

  -- ALGLIB --
     Copyright 03.09.2008 by Bochkanov Sergey
*************************************************************************/
double lrprocess(const linearmodel &lm, const real_1d_array &x);


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lrrmserror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lravgerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lravgrelerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);

/*************************************************************************
Filters: simple moving averages (unsymmetric).

This filter replaces array by results of SMA(K) filter. SMA(K) is defined
as filter which averages at most K previous points (previous - not points
AROUND central point) - or less, in case of the first K-1 points.

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    K           -   K>=1 (K can be larger than N ,  such  cases  will  be
                    correctly handled). Window width. K=1 corresponds  to
                    identity transformation (nothing changes).

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed with SMA(K)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm makes only one pass through array and uses running
        sum  to speed-up calculation of the averages. Additional measures
        are taken to ensure that running sum on a long sequence  of  zero
        elements will be correctly reset to zero even in the presence  of
        round-off error.

NOTE 3: this  is  unsymmetric version of the algorithm,  which  does  NOT
        averages points after the current one. Only X[i], X[i-1], ... are
        used when calculating new value of X[i]. We should also note that
        this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filtersma(real_1d_array &x, const ae_int_t n, const ae_int_t k);
void filtersma(real_1d_array &x, const ae_int_t k);


/*************************************************************************
Filters: exponential moving averages.

This filter replaces array by results of EMA(alpha) filter. EMA(alpha) is
defined as filter which replaces X[] by S[]:
    S[0] = X[0]
    S[t] = alpha*X[t] + (1-alpha)*S[t-1]

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    alpha       -   0<alpha<=1, smoothing parameter.

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed
                    with EMA(alpha)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

NOTE 3: technical analytis users quite often work  with  EMA  coefficient
        expressed in DAYS instead of fractions. If you want to  calculate
        EMA(N), where N is a number of days, you can use alpha=2/(N+1).

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filterema(real_1d_array &x, const ae_int_t n, const double alpha);
void filterema(real_1d_array &x, const double alpha);


/*************************************************************************
Filters: linear regression moving averages.

This filter replaces array by results of LRMA(K) filter.

LRMA(K) is defined as filter which, for each data  point,  builds  linear
regression  model  using  K  prevous  points (point itself is included in
these K points) and calculates value of this linear model at the point in
question.

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    K           -   K>=1 (K can be larger than N ,  such  cases  will  be
                    correctly handled). Window width. K=1 corresponds  to
                    identity transformation (nothing changes).

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed with SMA(K)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm makes only one pass through array and uses running
        sum  to speed-up calculation of the averages. Additional measures
        are taken to ensure that running sum on a long sequence  of  zero
        elements will be correctly reset to zero even in the presence  of
        round-off error.

NOTE 3: this  is  unsymmetric version of the algorithm,  which  does  NOT
        averages points after the current one. Only X[i], X[i-1], ... are
        used when calculating new value of X[i]. We should also note that
        this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filterlrma(real_1d_array &x, const ae_int_t n, const ae_int_t k);
void filterlrma(real_1d_array &x, const ae_int_t k);

/*************************************************************************
Multiclass Fisher LDA

Subroutine finds coefficients of linear combination which optimally separates
training set on classes.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two important  improvements   of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multithreading support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison. Best results are achieved  for  high-dimensional  problems
  ! (NVars is at least 256).
  !
  ! Multithreading is used to  accelerate  initial  phase  of  LDA,  which
  ! includes calculation of products of large matrices.  Again,  for  best
  ! efficiency problem must be high-dimensional.
  !
  ! Generally, commercial ALGLIB is several times faster than  open-source
  ! generic C edition, and many times faster than open-source C# edition.
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   linear combination coefficients, array[0..NVars-1]

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************/
void fisherlda(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_1d_array &w);


/*************************************************************************
N-dimensional multiclass Fisher LDA

Subroutine finds coefficients of linear combinations which optimally separates
training set on classes. It returns N-dimensional basis whose vector are sorted
by quality of training set separation (in descending order).

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes two important  improvements   of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  ! * multithreading support
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison. Best results are achieved  for  high-dimensional  problems
  ! (NVars is at least 256).
  !
  ! Multithreading is used to  accelerate  initial  phase  of  LDA,  which
  ! includes calculation of products of large matrices.  Again,  for  best
  ! efficiency problem must be high-dimensional.
  !
  ! Generally, commercial ALGLIB is several times faster than  open-source
  ! generic C edition, and many times faster than open-source C# edition.
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   basis, array[0..NVars-1,0..NVars-1]
                    columns of matrix stores basis vectors, sorted by
                    quality of training set separation (in descending order)

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************/
void fisherldan(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_2d_array &w);
void smp_fisherldan(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_2d_array &w);

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
void mlpserialize(multilayerperceptron &obj, std::string &s_out);


/*************************************************************************
This function unserializes data structure from string.
*************************************************************************/
void mlpunserialize(std::string &s_in, multilayerperceptron &obj);


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers, with linear output layer. Network weights are  filled  with  small
random values.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Same  as  MLPCreate0,  but  with  one  hidden  layer  (NHid  neurons) with
non-linear activation function. Output layer is linear.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreate0, but with two hidden layers (NHid1 and  NHid2  neurons)
with non-linear activation function. Output layer is linear.
 $ALL

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values.

Activation function of the output layer takes values:

    (B, +INF), if D>=0

or

    (-INF, B), if D<0.


  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateB0 but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateB0 but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values. Activation function of the output layer takes values [A,B].

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateR0, but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateR0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);


/*************************************************************************
Creates classifier network with NIn  inputs  and  NOut  possible  classes.
Network contains no hidden layers and linear output  layer  with  SOFTMAX-
normalization  (so  outputs  sums  up  to  1.0  and  converge to posterior
probabilities).

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateC0, but with one non-linear hidden layer.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateC0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Copying of neural network

INPUT PARAMETERS:
    Network1 -   original

OUTPUT PARAMETERS:
    Network2 -   copy

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcopy(const multilayerperceptron &network1, multilayerperceptron &network2);


/*************************************************************************
This function copies tunable  parameters (weights/means/sigmas)  from  one
network to another with same architecture. It  performs  some  rudimentary
checks that architectures are same, and throws exception if check fails.

It is intended for fast copying of states between two  network  which  are
known to have same geometry.

INPUT PARAMETERS:
    Network1 -   source, must be correctly initialized
    Network2 -   target, must have same architecture

OUTPUT PARAMETERS:
    Network2 -   network state is copied from source to target

  -- ALGLIB --
     Copyright 20.06.2013 by Bochkanov Sergey
*************************************************************************/
void mlpcopytunableparameters(const multilayerperceptron &network1, const multilayerperceptron &network2);


/*************************************************************************
Randomization of neural network weights

  -- ALGLIB --
     Copyright 06.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlprandomize(const multilayerperceptron &network);


/*************************************************************************
Randomization of neural network weights and standartisator

  -- ALGLIB --
     Copyright 10.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlprandomizefull(const multilayerperceptron &network);


/*************************************************************************
Internal subroutine.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpinitpreprocessor(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);


/*************************************************************************
Returns information about initialized network: number of inputs, outputs,
weights.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpproperties(const multilayerperceptron &network, ae_int_t &nin, ae_int_t &nout, ae_int_t &wcount);


/*************************************************************************
Returns number of inputs.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetinputscount(const multilayerperceptron &network);


/*************************************************************************
Returns number of outputs.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetoutputscount(const multilayerperceptron &network);


/*************************************************************************
Returns number of weights.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetweightscount(const multilayerperceptron &network);


/*************************************************************************
Tells whether network is SOFTMAX-normalized (i.e. classifier) or not.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
bool mlpissoftmax(const multilayerperceptron &network);


/*************************************************************************
This function returns total number of layers (including input, hidden and
output layers).

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetlayerscount(const multilayerperceptron &network);


/*************************************************************************
This function returns size of K-th layer.

K=0 corresponds to input layer, K=CNT-1 corresponds to output layer.

Size of the output layer is always equal to the number of outputs, although
when we have softmax-normalized network, last neuron doesn't have any
connections - it is just zero.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetlayersize(const multilayerperceptron &network, const ae_int_t k);


/*************************************************************************
This function returns offset/scaling coefficients for I-th input of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index

OUTPUT PARAMETERS:
    Mean        -   mean term
    Sigma       -   sigma term, guaranteed to be nonzero.

I-th input is passed through linear transformation
    IN[i] = (IN[i]-Mean)/Sigma
before feeding to the network

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpgetinputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma);


/*************************************************************************
This function returns offset/scaling coefficients for I-th output of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index

OUTPUT PARAMETERS:
    Mean        -   mean term
    Sigma       -   sigma term, guaranteed to be nonzero.

I-th output is passed through linear transformation
    OUT[i] = OUT[i]*Sigma+Mean
before returning it to user. In case we have SOFTMAX-normalized network,
we return (Mean,Sigma)=(0.0,1.0).

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpgetoutputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma);


/*************************************************************************
This function returns information about Ith neuron of Kth layer

INPUT PARAMETERS:
    Network     -   network
    K           -   layer index
    I           -   neuron index (within layer)

OUTPUT PARAMETERS:
    FKind       -   activation function type (used by MLPActivationFunction())
                    this value is zero for input or linear neurons
    Threshold   -   also called offset, bias
                    zero for input neurons

NOTE: this function throws exception if layer or neuron with  given  index
do not exists.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpgetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, ae_int_t &fkind, double &threshold);


/*************************************************************************
This function returns information about connection from I0-th neuron of
K0-th layer to I1-th neuron of K1-th layer.

INPUT PARAMETERS:
    Network     -   network
    K0          -   layer index
    I0          -   neuron index (within layer)
    K1          -   layer index
    I1          -   neuron index (within layer)

RESULT:
    connection weight (zero for non-existent connections)

This function:
1. throws exception if layer or neuron with given index do not exists.
2. returns zero if neurons exist, but there is no connection between them

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
double mlpgetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1);


/*************************************************************************
This function sets offset/scaling coefficients for I-th input of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index
    Mean        -   mean term
    Sigma       -   sigma term (if zero, will be replaced by 1.0)

NTE: I-th input is passed through linear transformation
    IN[i] = (IN[i]-Mean)/Sigma
before feeding to the network. This function sets Mean and Sigma.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetinputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma);


/*************************************************************************
This function sets offset/scaling coefficients for I-th output of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index
    Mean        -   mean term
    Sigma       -   sigma term (if zero, will be replaced by 1.0)

OUTPUT PARAMETERS:

NOTE: I-th output is passed through linear transformation
    OUT[i] = OUT[i]*Sigma+Mean
before returning it to user. This function sets Sigma/Mean. In case we
have SOFTMAX-normalized network, you can not set (Sigma,Mean) to anything
other than(0.0,1.0) - this function will throw exception.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetoutputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma);


/*************************************************************************
This function modifies information about Ith neuron of Kth layer

INPUT PARAMETERS:
    Network     -   network
    K           -   layer index
    I           -   neuron index (within layer)
    FKind       -   activation function type (used by MLPActivationFunction())
                    this value must be zero for input neurons
                    (you can not set activation function for input neurons)
    Threshold   -   also called offset, bias
                    this value must be zero for input neurons
                    (you can not set threshold for input neurons)

NOTES:
1. this function throws exception if layer or neuron with given index do
   not exists.
2. this function also throws exception when you try to set non-linear
   activation function for input neurons (any kind of network) or for output
   neurons of classifier network.
3. this function throws exception when you try to set non-zero threshold for
   input neurons (any kind of network).

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, const ae_int_t fkind, const double threshold);


/*************************************************************************
This function modifies information about connection from I0-th neuron of
K0-th layer to I1-th neuron of K1-th layer.

INPUT PARAMETERS:
    Network     -   network
    K0          -   layer index
    I0          -   neuron index (within layer)
    K1          -   layer index
    I1          -   neuron index (within layer)
    W           -   connection weight (must be zero for non-existent
                    connections)

This function:
1. throws exception if layer or neuron with given index do not exists.
2. throws exception if you try to set non-zero weight for non-existent
   connection

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1, const double w);


/*************************************************************************
Neural network activation function

INPUT PARAMETERS:
    NET         -   neuron input
    K           -   function index (zero for linear function)

OUTPUT PARAMETERS:
    F           -   function
    DF          -   its derivative
    D2F         -   its second derivative

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpactivationfunction(const double net, const ae_int_t k, double &f, double &df, double &d2f);


/*************************************************************************
Procesing

INPUT PARAMETERS:
    Network -   neural network
    X       -   input vector,  array[0..NIn-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

See also MLPProcessI

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpprocess(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y);


/*************************************************************************
'interactive'  variant  of  MLPProcess  for  languages  like  Python which
support constructs like "Y = MLPProcess(NN,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 21.09.2010 by Bochkanov Sergey
*************************************************************************/
void mlpprocessi(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y);


/*************************************************************************
Error of the neural network on dataset.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore systems.
  ! Second improvement gives constant speedup (2-3x, depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format;
    NPoints     -   points count.

RESULT:
    sum-of-squares error, SUM(sqr(y[i]-desired_y[i])/2)

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double smp_mlperror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Error of the neural network on dataset given by sparse matrix.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore systems.
  ! Second improvement gives constant speedup (2-3x, depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network
    XY          -   training  set,  see  below  for  information  on   the
                    training set format. This function checks  correctness
                    of  the  dataset  (no  NANs/INFs,  class  numbers  are
                    correct) and throws exception when  incorrect  dataset
                    is passed.  Sparse  matrix  must  use  CRS  format for
                    storage.
    NPoints     -   points count, >=0

RESULT:
    sum-of-squares error, SUM(sqr(y[i]-desired_y[i])/2)

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
double mlperrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double smp_mlperrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);


/*************************************************************************
Natural error function for neural network, internal subroutine.

NOTE: this function is single-threaded. Unlike other  error  function,  it
receives no speed-up from being executed in SMP mode.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperrorn(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);


/*************************************************************************
Classification error of the neural network on dataset.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format;
    NPoints     -   points count.

RESULT:
    classification error (number of misclassified cases)

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
ae_int_t smp_mlpclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Relative classification error on the test set.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format;
    NPoints     -   points count.

RESULT:
Percent   of incorrectly   classified  cases.  Works  both  for classifier
networks and general purpose networks used as classifiers.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 25.12.2008 by Bochkanov Sergey
*************************************************************************/
double mlprelclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double smp_mlprelclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Relative classification error on the test set given by sparse matrix.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format. Sparse matrix must use CRS format
                    for storage.
    NPoints     -   points count, >=0.

RESULT:
Percent   of incorrectly   classified  cases.  Works  both  for classifier
networks and general purpose networks used as classifiers.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 09.08.2012 by Bochkanov Sergey
*************************************************************************/
double mlprelclserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double smp_mlprelclserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);


/*************************************************************************
Average cross-entropy  (in bits  per element) on the test set.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format;
    NPoints     -   points count.

RESULT:
CrossEntropy/(NPoints*LN(2)).
Zero if network solves regression task.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 08.01.2009 by Bochkanov Sergey
*************************************************************************/
double mlpavgce(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double smp_mlpavgce(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average  cross-entropy  (in bits  per element)  on the  test set  given by
sparse matrix.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format. This function checks  correctness
                    of  the  dataset  (no  NANs/INFs,  class  numbers  are
                    correct) and throws exception when  incorrect  dataset
                    is passed.  Sparse  matrix  must  use  CRS  format for
                    storage.
    NPoints     -   points count, >=0.

RESULT:
CrossEntropy/(NPoints*LN(2)).
Zero if network solves regression task.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 9.08.2012 by Bochkanov Sergey
*************************************************************************/
double mlpavgcesparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double smp_mlpavgcesparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);


/*************************************************************************
RMS error on the test set given.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format;
    NPoints     -   points count.

RESULT:
Root mean  square error. Its meaning for regression task is obvious. As for
classification  task,  RMS  error  means  error  when estimating  posterior
probabilities.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlprmserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double smp_mlprmserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
RMS error on the test set given by sparse matrix.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format. This function checks  correctness
                    of  the  dataset  (no  NANs/INFs,  class  numbers  are
                    correct) and throws exception when  incorrect  dataset
                    is passed.  Sparse  matrix  must  use  CRS  format for
                    storage.
    NPoints     -   points count, >=0.

RESULT:
Root mean  square error. Its meaning for regression task is obvious. As for
classification  task,  RMS  error  means  error  when estimating  posterior
probabilities.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 09.08.2012 by Bochkanov Sergey
*************************************************************************/
double mlprmserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double smp_mlprmserrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);


/*************************************************************************
Average absolute error on the test set.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format;
    NPoints     -   points count.

RESULT:
Its meaning for regression task is obvious. As for classification task, it
means average error when estimating posterior probabilities.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 11.03.2008 by Bochkanov Sergey
*************************************************************************/
double mlpavgerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double smp_mlpavgerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average absolute error on the test set given by sparse matrix.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format. This function checks  correctness
                    of  the  dataset  (no  NANs/INFs,  class  numbers  are
                    correct) and throws exception when  incorrect  dataset
                    is passed.  Sparse  matrix  must  use  CRS  format for
                    storage.
    NPoints     -   points count, >=0.

RESULT:
Its meaning for regression task is obvious. As for classification task, it
means average error when estimating posterior probabilities.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 09.08.2012 by Bochkanov Sergey
*************************************************************************/
double mlpavgerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double smp_mlpavgerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);


/*************************************************************************
Average relative error on the test set.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format;
    NPoints     -   points count.

RESULT:
Its meaning for regression task is obvious. As for classification task, it
means  average  relative  error  when  estimating posterior probability of
belonging to the correct class.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 11.03.2008 by Bochkanov Sergey
*************************************************************************/
double mlpavgrelerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);
double smp_mlpavgrelerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average relative error on the test set given by sparse matrix.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network     -   neural network;
    XY          -   training  set,  see  below  for  information  on   the
                    training set format. This function checks  correctness
                    of  the  dataset  (no  NANs/INFs,  class  numbers  are
                    correct) and throws exception when  incorrect  dataset
                    is passed.  Sparse  matrix  must  use  CRS  format for
                    storage.
    NPoints     -   points count, >=0.

RESULT:
Its meaning for regression task is obvious. As for classification task, it
means  average  relative  error  when  estimating posterior probability of
belonging to the correct class.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 09.08.2012 by Bochkanov Sergey
*************************************************************************/
double mlpavgrelerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);
double smp_mlpavgrelerrorsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t npoints);


/*************************************************************************
Gradient calculation

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    X       -   input vector, length of array must be at least NIn
    DesiredY-   desired outputs, length of array must be at least NOut
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgrad(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad);


/*************************************************************************
Gradient calculation (natural error function is used)

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    X       -   input vector, length of array must be at least NIn
    DesiredY-   desired outputs, length of array must be at least NOut
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, sum-of-squares for regression networks,
                cross-entropy for classification networks.
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradn(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad);


/*************************************************************************
Batch gradient calculation for a set of inputs/outputs


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   original dataset in dense format; one sample = one row:
                * first NIn columns contain inputs,
                * for regression problem, next NOut columns store
                  desired outputs.
                * for classification problem, next column (just one!)
                  stores class number.
    SSize   -   number of elements in XY
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);
void smp_mlpgradbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);


/*************************************************************************
Batch gradient calculation for a set  of inputs/outputs  given  by  sparse
matrices


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   original dataset in sparse format; one sample = one row:
                * MATRIX MUST BE STORED IN CRS FORMAT
                * first NIn columns contain inputs.
                * for regression problem, next NOut columns store
                  desired outputs.
                * for classification problem, next column (just one!)
                  stores class number.
    SSize   -   number of elements in XY
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 26.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpgradbatchsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t ssize, double &e, real_1d_array &grad);
void smp_mlpgradbatchsparse(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t ssize, double &e, real_1d_array &grad);


/*************************************************************************
Batch gradient calculation for a subset of dataset


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   original dataset in dense format; one sample = one row:
                * first NIn columns contain inputs,
                * for regression problem, next NOut columns store
                  desired outputs.
                * for classification problem, next column (just one!)
                  stores class number.
    SetSize -   real size of XY, SetSize>=0;
    Idx     -   subset of SubsetSize elements, array[SubsetSize]:
                * Idx[I] stores row index in the original dataset which is
                  given by XY. Gradient is calculated with respect to rows
                  whose indexes are stored in Idx[].
                * Idx[]  must store correct indexes; this function  throws
                  an  exception  in  case  incorrect index (less than 0 or
                  larger than rows(XY)) is given
                * Idx[]  may  store  indexes  in  any  order and even with
                  repetitions.
    SubsetSize- number of elements in Idx[] array:
                * positive value means that subset given by Idx[] is processed
                * zero value results in zero gradient
                * negative value means that full dataset is processed
    Grad      - possibly  preallocated array. If size of array is  smaller
                than WCount, it will be reallocated. It is  recommended to
                reuse  previously  allocated  array  to  reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E         - error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad      - gradient  of  E  with  respect   to  weights  of  network,
                array[WCount]

  -- ALGLIB --
     Copyright 26.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpgradbatchsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);
void smp_mlpgradbatchsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);


/*************************************************************************
Batch gradient calculation for a set of inputs/outputs  for  a  subset  of
dataset given by set of indexes.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   original dataset in sparse format; one sample = one row:
                * MATRIX MUST BE STORED IN CRS FORMAT
                * first NIn columns contain inputs,
                * for regression problem, next NOut columns store
                  desired outputs.
                * for classification problem, next column (just one!)
                  stores class number.
    SetSize -   real size of XY, SetSize>=0;
    Idx     -   subset of SubsetSize elements, array[SubsetSize]:
                * Idx[I] stores row index in the original dataset which is
                  given by XY. Gradient is calculated with respect to rows
                  whose indexes are stored in Idx[].
                * Idx[]  must store correct indexes; this function  throws
                  an  exception  in  case  incorrect index (less than 0 or
                  larger than rows(XY)) is given
                * Idx[]  may  store  indexes  in  any  order and even with
                  repetitions.
    SubsetSize- number of elements in Idx[] array:
                * positive value means that subset given by Idx[] is processed
                * zero value results in zero gradient
                * negative value means that full dataset is processed
    Grad      - possibly  preallocated array. If size of array is  smaller
                than WCount, it will be reallocated. It is  recommended to
                reuse  previously  allocated  array  to  reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad    -   gradient  of  E  with  respect   to  weights  of  network,
                array[WCount]

NOTE: when  SubsetSize<0 is used full dataset by call MLPGradBatchSparse
      function.

  -- ALGLIB --
     Copyright 26.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpgradbatchsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);
void smp_mlpgradbatchsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &idx, const ae_int_t subsetsize, double &e, real_1d_array &grad);


/*************************************************************************
Batch gradient calculation for a set of inputs/outputs
(natural error function is used)

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   set of inputs/outputs; one sample = one row;
                first NIn columns contain inputs,
                next NOut columns - desired outputs.
    SSize   -   number of elements in XY
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, sum-of-squares for regression networks,
                cross-entropy for classification networks.
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradnbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);


/*************************************************************************
Batch Hessian calculation (natural error function) using R-algorithm.
Internal subroutine.

  -- ALGLIB --
     Copyright 26.01.2008 by Bochkanov Sergey.

     Hessian calculation based on R-algorithm described in
     "Fast Exact Multiplication by the Hessian",
     B. A. Pearlmutter,
     Neural Computation, 1994.
*************************************************************************/
void mlphessiannbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h);


/*************************************************************************
Batch Hessian calculation using R-algorithm.
Internal subroutine.

  -- ALGLIB --
     Copyright 26.01.2008 by Bochkanov Sergey.

     Hessian calculation based on R-algorithm described in
     "Fast Exact Multiplication by the Hessian",
     B. A. Pearlmutter,
     Neural Computation, 1994.
*************************************************************************/
void mlphessianbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h);


/*************************************************************************
Calculation of all types of errors on subset of dataset.

FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   original dataset; one sample = one row;
                first NIn columns contain inputs,
                next NOut columns - desired outputs.
    SetSize -   real size of XY, SetSize>=0;
    Subset  -   subset of SubsetSize elements, array[SubsetSize];
    SubsetSize- number of elements in Subset[] array:
                * if SubsetSize>0, rows of XY with indices Subset[0]...
                  ...Subset[SubsetSize-1] are processed
                * if SubsetSize=0, zeros are returned
                * if SubsetSize<0, entire dataset is  processed;  Subset[]
                  array is ignored in this case.

OUTPUT PARAMETERS:
    Rep     -   it contains all type of errors.

  -- ALGLIB --
     Copyright 04.09.2012 by Bochkanov Sergey
*************************************************************************/
void mlpallerrorssubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);
void smp_mlpallerrorssubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);


/*************************************************************************
Calculation of all types of errors on subset of dataset.

FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   original dataset given by sparse matrix;
                one sample = one row;
                first NIn columns contain inputs,
                next NOut columns - desired outputs.
    SetSize -   real size of XY, SetSize>=0;
    Subset  -   subset of SubsetSize elements, array[SubsetSize];
    SubsetSize- number of elements in Subset[] array:
                * if SubsetSize>0, rows of XY with indices Subset[0]...
                  ...Subset[SubsetSize-1] are processed
                * if SubsetSize=0, zeros are returned
                * if SubsetSize<0, entire dataset is  processed;  Subset[]
                  array is ignored in this case.

OUTPUT PARAMETERS:
    Rep     -   it contains all type of errors.


  -- ALGLIB --
     Copyright 04.09.2012 by Bochkanov Sergey
*************************************************************************/
void mlpallerrorssparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);
void smp_mlpallerrorssparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize, modelerrors &rep);


/*************************************************************************
Error of the neural network on subset of dataset.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network   -     neural network;
    XY        -     training  set,  see  below  for  information  on   the
                    training set format;
    SetSize   -     real size of XY, SetSize>=0;
    Subset    -     subset of SubsetSize elements, array[SubsetSize];
    SubsetSize-     number of elements in Subset[] array:
                    * if SubsetSize>0, rows of XY with indices Subset[0]...
                      ...Subset[SubsetSize-1] are processed
                    * if SubsetSize=0, zeros are returned
                    * if SubsetSize<0, entire dataset is  processed;  Subset[]
                      array is ignored in this case.

RESULT:
    sum-of-squares error, SUM(sqr(y[i]-desired_y[i])/2)

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 04.09.2012 by Bochkanov Sergey
*************************************************************************/
double mlperrorsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);
double smp_mlperrorsubset(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);


/*************************************************************************
Error of the neural network on subset of sparse dataset.


FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support
  !
  ! First improvement gives close-to-linear speedup on multicore  systems.
  ! Second improvement gives constant speedup (2-3x depending on your CPU)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.


INPUT PARAMETERS:
    Network   -     neural network;
    XY        -     training  set,  see  below  for  information  on   the
                    training set format. This function checks  correctness
                    of  the  dataset  (no  NANs/INFs,  class  numbers  are
                    correct) and throws exception when  incorrect  dataset
                    is passed.  Sparse  matrix  must  use  CRS  format for
                    storage.
    SetSize   -     real size of XY, SetSize>=0;
                    it is used when SubsetSize<0;
    Subset    -     subset of SubsetSize elements, array[SubsetSize];
    SubsetSize-     number of elements in Subset[] array:
                    * if SubsetSize>0, rows of XY with indices Subset[0]...
                      ...Subset[SubsetSize-1] are processed
                    * if SubsetSize=0, zeros are returned
                    * if SubsetSize<0, entire dataset is  processed;  Subset[]
                      array is ignored in this case.

RESULT:
    sum-of-squares error, SUM(sqr(y[i]-desired_y[i])/2)

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
dataset format is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 04.09.2012 by Bochkanov Sergey
*************************************************************************/
double mlperrorsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);
double smp_mlperrorsparsesubset(const multilayerperceptron &network, const sparsematrix &xy, const ae_int_t setsize, const integer_1d_array &subset, const ae_int_t subsetsize);

/*************************************************************************
This subroutine trains logit model.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars]
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<NVars+2, NVars<1, NClasses<2).
                    *  1, if task has been solved
    LM          -   model built
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnltrainh(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, logitmodel &lm, mnlreport &rep);


/*************************************************************************
Procesing

INPUT PARAMETERS:
    LM      -   logit model, passed by non-constant reference
                (some fields of structure are used as temporaries
                when calculating model output).
    X       -   input vector,  array[0..NVars-1].
    Y       -   (possibly) preallocated buffer; if size of Y is less than
                NClasses, it will be reallocated.If it is large enough, it
                is NOT reallocated, so we can save some time on reallocation.

OUTPUT PARAMETERS:
    Y       -   result, array[0..NClasses-1]
                Vector of posterior probabilities for classification task.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlprocess(const logitmodel &lm, const real_1d_array &x, real_1d_array &y);


/*************************************************************************
'interactive'  variant  of  MNLProcess  for  languages  like  Python which
support constructs like "Y = MNLProcess(LM,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlprocessi(const logitmodel &lm, const real_1d_array &x, real_1d_array &y);


/*************************************************************************
Unpacks coefficients of logit model. Logit model have form:

    P(class=i) = S(i) / (S(0) + S(1) + ... +S(M-1))
          S(i) = Exp(A[i,0]*X[0] + ... + A[i,N-1]*X[N-1] + A[i,N]), when i<M-1
        S(M-1) = 1

INPUT PARAMETERS:
    LM          -   logit model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NClasses-2,0..NVars]
    NVars       -   number of independent variables
    NClasses    -   number of classes

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlunpack(const logitmodel &lm, real_2d_array &a, ae_int_t &nvars, ae_int_t &nclasses);


/*************************************************************************
"Packs" coefficients and creates logit model in ALGLIB format (MNLUnpack
reversed).

INPUT PARAMETERS:
    A           -   model (see MNLUnpack)
    NVars       -   number of independent variables
    NClasses    -   number of classes

OUTPUT PARAMETERS:
    LM          -   logit model.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlpack(const real_2d_array &a, const ae_int_t nvars, const ae_int_t nclasses, logitmodel &lm);


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*ln(2)).

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgce(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
double mnlrelclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlrmserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgrelerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t ssize);


/*************************************************************************
Classification error on test set = MNLRelClsError*NPoints

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
ae_int_t mnlclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);

/*************************************************************************
DESCRIPTION:

This function creates MCPD (Markov Chains for Population Data) solver.

This  solver  can  be  used  to find transition matrix P for N-dimensional
prediction  problem  where transition from X[i] to X[i+1] is  modelled  as
    X[i+1] = P*X[i]
where X[i] and X[i+1] are N-dimensional population vectors (components  of
each X are non-negative), and P is a N*N transition matrix (elements of  P
are non-negative, each column sums to 1.0).

Such models arise when when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is constant, i.e. there is no new individuals and no one
  leaves population
* you want to model transitions of individuals from one state into another

USAGE:

Here we give very brief outline of the MCPD. We strongly recommend you  to
read examples in the ALGLIB Reference Manual and to read ALGLIB User Guide
on data analysis which is available at http://www.alglib.net/dataanalysis/

1. User initializes algorithm state with MCPDCreate() call

2. User  adds  one  or  more  tracks -  sequences of states which describe
   evolution of a system being modelled from different starting conditions

3. User may add optional boundary, equality  and/or  linear constraints on
   the coefficients of P by calling one of the following functions:
   * MCPDSetEC() to set equality constraints
   * MCPDSetBC() to set bound constraints
   * MCPDSetLC() to set linear constraints

4. Optionally,  user  may  set  custom  weights  for prediction errors (by
   default, algorithm assigns non-equal, automatically chosen weights  for
   errors in the prediction of different components of X). It can be  done
   with a call of MCPDSetPredictionWeights() function.

5. User calls MCPDSolve() function which takes algorithm  state and
   pointer (delegate, etc.) to callback function which calculates F/G.

6. User calls MCPDResults() to get solution

INPUT PARAMETERS:
    N       -   problem dimension, N>=1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreate(const ae_int_t n, mcpdstate &s);


/*************************************************************************
DESCRIPTION:

This function is a specialized version of MCPDCreate()  function,  and  we
recommend  you  to read comments for this function for general information
about MCPD solver.

This  function  creates  MCPD (Markov Chains for Population  Data)  solver
for "Entry-state" model,  i.e. model  where transition from X[i] to X[i+1]
is modelled as
    X[i+1] = P*X[i]
where
    X[i] and X[i+1] are N-dimensional state vectors
    P is a N*N transition matrix
and  one  selected component of X[] is called "entry" state and is treated
in a special way:
    system state always transits from "entry" state to some another state
    system state can not transit from any state into "entry" state
Such conditions basically mean that row of P which corresponds to  "entry"
state is zero.

Such models arise when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is NOT constant -  at every moment of time there is some
  (unpredictable) amount of "new" individuals, which can transit into  one
  of the states at the next turn, but still no one leaves population
* you want to model transitions of individuals from one state into another
* but you do NOT want to predict amount of "new"  individuals  because  it
  does not depends on individuals already present (hence  system  can  not
  transit INTO entry state - it can only transit FROM it).

This model is discussed  in  more  details  in  the ALGLIB User Guide (see
http://www.alglib.net/dataanalysis/ for more data).

INPUT PARAMETERS:
    N       -   problem dimension, N>=2
    EntryState- index of entry state, in 0..N-1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreateentry(const ae_int_t n, const ae_int_t entrystate, mcpdstate &s);


/*************************************************************************
DESCRIPTION:

This function is a specialized version of MCPDCreate()  function,  and  we
recommend  you  to read comments for this function for general information
about MCPD solver.

This  function  creates  MCPD (Markov Chains for Population  Data)  solver
for "Exit-state" model,  i.e. model  where  transition from X[i] to X[i+1]
is modelled as
    X[i+1] = P*X[i]
where
    X[i] and X[i+1] are N-dimensional state vectors
    P is a N*N transition matrix
and  one  selected component of X[] is called "exit"  state and is treated
in a special way:
    system state can transit from any state into "exit" state
    system state can not transit from "exit" state into any other state
    transition operator discards "exit" state (makes it zero at each turn)
Such  conditions  basically  mean  that  column  of P which corresponds to
"exit" state is zero. Multiplication by such P may decrease sum of  vector
components.

Such models arise when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is NOT constant - individuals can move into "exit" state
  and leave population at the next turn, but there are no new individuals
* amount of individuals which leave population can be predicted
* you want to model transitions of individuals from one state into another
  (including transitions into the "exit" state)

This model is discussed  in  more  details  in  the ALGLIB User Guide (see
http://www.alglib.net/dataanalysis/ for more data).

INPUT PARAMETERS:
    N       -   problem dimension, N>=2
    ExitState-  index of exit state, in 0..N-1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreateexit(const ae_int_t n, const ae_int_t exitstate, mcpdstate &s);


/*************************************************************************
DESCRIPTION:

This function is a specialized version of MCPDCreate()  function,  and  we
recommend  you  to read comments for this function for general information
about MCPD solver.

This  function  creates  MCPD (Markov Chains for Population  Data)  solver
for "Entry-Exit-states" model, i.e. model where  transition  from  X[i] to
X[i+1] is modelled as
    X[i+1] = P*X[i]
where
    X[i] and X[i+1] are N-dimensional state vectors
    P is a N*N transition matrix
one selected component of X[] is called "entry" state and is treated in  a
special way:
    system state always transits from "entry" state to some another state
    system state can not transit from any state into "entry" state
and another one component of X[] is called "exit" state and is treated  in
a special way too:
    system state can transit from any state into "exit" state
    system state can not transit from "exit" state into any other state
    transition operator discards "exit" state (makes it zero at each turn)
Such conditions basically mean that:
    row of P which corresponds to "entry" state is zero
    column of P which corresponds to "exit" state is zero
Multiplication by such P may decrease sum of vector components.

Such models arise when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is NOT constant
* at every moment of time there is some (unpredictable)  amount  of  "new"
  individuals, which can transit into one of the states at the next turn
* some  individuals  can  move  (predictably)  into "exit" state and leave
  population at the next turn
* you want to model transitions of individuals from one state into another,
  including transitions from the "entry" state and into the "exit" state.
* but you do NOT want to predict amount of "new"  individuals  because  it
  does not depends on individuals already present (hence  system  can  not
  transit INTO entry state - it can only transit FROM it).

This model is discussed  in  more  details  in  the ALGLIB User Guide (see
http://www.alglib.net/dataanalysis/ for more data).

INPUT PARAMETERS:
    N       -   problem dimension, N>=2
    EntryState- index of entry state, in 0..N-1
    ExitState-  index of exit state, in 0..N-1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreateentryexit(const ae_int_t n, const ae_int_t entrystate, const ae_int_t exitstate, mcpdstate &s);


/*************************************************************************
This  function  is  used to add a track - sequence of system states at the
different moments of its evolution.

You  may  add  one  or several tracks to the MCPD solver. In case you have
several tracks, they won't overwrite each other. For example,  if you pass
two tracks, A1-A2-A3 (system at t=A+1, t=A+2 and t=A+3) and B1-B2-B3, then
solver will try to model transitions from t=A+1 to t=A+2, t=A+2 to  t=A+3,
t=B+1 to t=B+2, t=B+2 to t=B+3. But it WONT mix these two tracks - i.e. it
wont try to model transition from t=A+3 to t=B+1.

INPUT PARAMETERS:
    S       -   solver
    XY      -   track, array[K,N]:
                * I-th row is a state at t=I
                * elements of XY must be non-negative (exception will be
                  thrown on negative elements)
    K       -   number of points in a track
                * if given, only leading K rows of XY are used
                * if not given, automatically determined from size of XY

NOTES:

1. Track may contain either proportional or population data:
   * with proportional data all rows of XY must sum to 1.0, i.e. we have
     proportions instead of absolute population values
   * with population data rows of XY contain population counts and generally
     do not sum to 1.0 (although they still must be non-negative)

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy, const ae_int_t k);
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy);


/*************************************************************************
This function is used to add equality constraints on the elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This  function  can  be  used  to  place equality constraints on arbitrary
subset of elements of P. Set of constraints is specified by EC, which  may
contain either NAN's or finite numbers from [0,1]. NAN denotes absence  of
constraint, finite number denotes equality constraint on specific  element
of P.

You can also  use  MCPDAddEC()  function  which  allows  to  ADD  equality
constraint  for  one  element  of P without changing constraints for other
elements.

These functions (MCPDSetEC and MCPDAddEC) interact as follows:
* there is internal matrix of equality constraints which is stored in  the
  MCPD solver
* MCPDSetEC() replaces this matrix by another one (SET)
* MCPDAddEC() modifies one element of this matrix and  leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddEC()  call  preserves  all  modifications  done by previous
  calls,  while  MCPDSetEC()  completely discards all changes  done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    EC      -   equality constraints, array[N,N]. Elements of  EC  can  be
                either NAN's or finite  numbers from  [0,1].  NAN  denotes
                absence  of  constraints,  while  finite  value    denotes
                equality constraint on the corresponding element of P.

NOTES:

1. infinite values of EC will lead to exception being thrown. Values  less
than 0.0 or greater than 1.0 will lead to error code being returned  after
call to MCPDSolve().

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetec(const mcpdstate &s, const real_2d_array &ec);


/*************************************************************************
This function is used to add equality constraints on the elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This function can be used to ADD equality constraint for one element of  P
without changing constraints for other elements.

You  can  also  use  MCPDSetEC()  function  which  allows  you  to specify
arbitrary set of equality constraints in one call.

These functions (MCPDSetEC and MCPDAddEC) interact as follows:
* there is internal matrix of equality constraints which is stored in the
  MCPD solver
* MCPDSetEC() replaces this matrix by another one (SET)
* MCPDAddEC() modifies one element of this matrix and leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddEC()  call  preserves  all  modifications done by previous
  calls,  while  MCPDSetEC()  completely discards all changes done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    I       -   row index of element being constrained
    J       -   column index of element being constrained
    C       -   value (constraint for P[I,J]).  Can  be  either  NAN  (no
                constraint) or finite value from [0,1].

NOTES:

1. infinite values of C  will lead to exception being thrown. Values  less
than 0.0 or greater than 1.0 will lead to error code being returned  after
call to MCPDSolve().

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdaddec(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double c);


/*************************************************************************
This function is used to add bound constraints  on  the  elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This  function  can  be  used  to  place bound   constraints  on arbitrary
subset  of  elements  of  P.  Set of constraints is specified by BndL/BndU
matrices, which may contain arbitrary combination  of  finite  numbers  or
infinities (like -INF<x<=0.5 or 0.1<=x<+INF).

You can also use MCPDAddBC() function which allows to ADD bound constraint
for one element of P without changing constraints for other elements.

These functions (MCPDSetBC and MCPDAddBC) interact as follows:
* there is internal matrix of bound constraints which is stored in the
  MCPD solver
* MCPDSetBC() replaces this matrix by another one (SET)
* MCPDAddBC() modifies one element of this matrix and  leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddBC()  call  preserves  all  modifications  done by previous
  calls,  while  MCPDSetBC()  completely discards all changes  done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    BndL    -   lower bounds constraints, array[N,N]. Elements of BndL can
                be finite numbers or -INF.
    BndU    -   upper bounds constraints, array[N,N]. Elements of BndU can
                be finite numbers or +INF.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetbc(const mcpdstate &s, const real_2d_array &bndl, const real_2d_array &bndu);


/*************************************************************************
This function is used to add bound constraints  on  the  elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This  function  can  be  used to ADD bound constraint for one element of P
without changing constraints for other elements.

You  can  also  use  MCPDSetBC()  function  which  allows to  place  bound
constraints  on arbitrary subset of elements of P.   Set of constraints is
specified  by  BndL/BndU matrices, which may contain arbitrary combination
of finite numbers or infinities (like -INF<x<=0.5 or 0.1<=x<+INF).

These functions (MCPDSetBC and MCPDAddBC) interact as follows:
* there is internal matrix of bound constraints which is stored in the
  MCPD solver
* MCPDSetBC() replaces this matrix by another one (SET)
* MCPDAddBC() modifies one element of this matrix and  leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddBC()  call  preserves  all  modifications  done by previous
  calls,  while  MCPDSetBC()  completely discards all changes  done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    I       -   row index of element being constrained
    J       -   column index of element being constrained
    BndL    -   lower bound
    BndU    -   upper bound

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdaddbc(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double bndl, const double bndu);


/*************************************************************************
This function is used to set linear equality/inequality constraints on the
elements of the transition matrix P.

This function can be used to set one or several general linear constraints
on the elements of P. Two types of constraints are supported:
* equality constraints
* inequality constraints (both less-or-equal and greater-or-equal)

Coefficients  of  constraints  are  specified  by  matrix  C (one  of  the
parameters).  One  row  of  C  corresponds  to  one  constraint.   Because
transition  matrix P has N*N elements,  we  need  N*N columns to store all
coefficients  (they  are  stored row by row), and one more column to store
right part - hence C has N*N+1 columns.  Constraint  kind is stored in the
CT array.

Thus, I-th linear constraint is
    P[0,0]*C[I,0] + P[0,1]*C[I,1] + .. + P[0,N-1]*C[I,N-1] +
        + P[1,0]*C[I,N] + P[1,1]*C[I,N+1] + ... +
        + P[N-1,N-1]*C[I,N*N-1]  ?=?  C[I,N*N]
where ?=? can be either "=" (CT[i]=0), "<=" (CT[i]<0) or ">=" (CT[i]>0).

Your constraint may involve only some subset of P (less than N*N elements).
For example it can be something like
    P[0,0] + P[0,1] = 0.5
In this case you still should pass matrix  with N*N+1 columns, but all its
elements (except for C[0,0], C[0,1] and C[0,N*N-1]) will be zero.

INPUT PARAMETERS:
    S       -   solver
    C       -   array[K,N*N+1] - coefficients of constraints
                (see above for complete description)
    CT      -   array[K] - constraint types
                (see above for complete description)
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k);
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct);


/*************************************************************************
This function allows to  tune  amount  of  Tikhonov  regularization  being
applied to your problem.

By default, regularizing term is equal to r*||P-prior_P||^2, where r is  a
small non-zero value,  P is transition matrix, prior_P is identity matrix,
||X||^2 is a sum of squared elements of X.

This  function  allows  you to change coefficient r. You can  also  change
prior values with MCPDSetPrior() function.

INPUT PARAMETERS:
    S       -   solver
    V       -   regularization  coefficient, finite non-negative value. It
                is  not  recommended  to specify zero value unless you are
                pretty sure that you want it.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsettikhonovregularizer(const mcpdstate &s, const double v);


/*************************************************************************
This  function  allows to set prior values used for regularization of your
problem.

By default, regularizing term is equal to r*||P-prior_P||^2, where r is  a
small non-zero value,  P is transition matrix, prior_P is identity matrix,
||X||^2 is a sum of squared elements of X.

This  function  allows  you to change prior values prior_P. You  can  also
change r with MCPDSetTikhonovRegularizer() function.

INPUT PARAMETERS:
    S       -   solver
    PP      -   array[N,N], matrix of prior values:
                1. elements must be real numbers from [0,1]
                2. columns must sum to 1.0.
                First property is checked (exception is thrown otherwise),
                while second one is not checked/enforced.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetprior(const mcpdstate &s, const real_2d_array &pp);


/*************************************************************************
This function is used to change prediction weights

MCPD solver scales prediction errors as follows
    Error(P) = ||W*(y-P*x)||^2
where
    x is a system state at time t
    y is a system state at time t+1
    P is a transition matrix
    W is a diagonal scaling matrix

By default, weights are chosen in order  to  minimize  relative prediction
error instead of absolute one. For example, if one component of  state  is
about 0.5 in magnitude and another one is about 0.05, then algorithm  will
make corresponding weights equal to 2.0 and 20.0.

INPUT PARAMETERS:
    S       -   solver
    PW      -   array[N], weights:
                * must be non-negative values (exception will be thrown otherwise)
                * zero values will be replaced by automatically chosen values

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetpredictionweights(const mcpdstate &s, const real_1d_array &pw);


/*************************************************************************
This function is used to start solution of the MCPD problem.

After return from this function, you can use MCPDResults() to get solution
and completion code.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsolve(const mcpdstate &s);


/*************************************************************************
MCPD results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    P       -   array[N,N], transition matrix
    Rep     -   optimization report. You should check Rep.TerminationType
                in  order  to  distinguish  successful  termination  from
                unsuccessful one. Speaking short, positive values  denote
                success, negative ones are failures.
                More information about fields of this  structure  can  be
                found in the comments on MCPDReport datatype.


  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdresults(const mcpdstate &s, real_2d_array &p, mcpdreport &rep);

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
void mlpeserialize(mlpensemble &obj, std::string &s_out);


/*************************************************************************
This function unserializes data structure from string.
*************************************************************************/
void mlpeunserialize(std::string &s_in, mlpensemble &obj);


/*************************************************************************
Like MLPCreate0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreate1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreate2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateB0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateB1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateB2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateR0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateR1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateR2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateC0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateC1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateC2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Creates ensemble from network. Only network geometry is copied.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatefromnetwork(const multilayerperceptron &network, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Randomization of MLP ensemble

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlperandomize(const mlpensemble &ensemble);


/*************************************************************************
Return ensemble properties (number of inputs and outputs).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeproperties(const mlpensemble &ensemble, ae_int_t &nin, ae_int_t &nout);


/*************************************************************************
Return normalization type (whether ensemble is SOFTMAX-normalized or not).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
bool mlpeissoftmax(const mlpensemble &ensemble);


/*************************************************************************
Procesing

INPUT PARAMETERS:
    Ensemble-   neural networks ensemble
    X       -   input vector,  array[0..NIn-1].
    Y       -   (possibly) preallocated buffer; if size of Y is less than
                NOut, it will be reallocated. If it is large enough, it
                is NOT reallocated, so we can save some time on reallocation.


OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeprocess(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y);


/*************************************************************************
'interactive'  variant  of  MLPEProcess  for  languages  like Python which
support constructs like "Y = MLPEProcess(LM,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeprocessi(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y);


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Works both for classifier betwork and for regression networks which
are used as classifiers.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlperelclserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if ensemble solves regression task.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgce(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for classification task
RMS error means error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpermserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average relative error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgrelerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);

/*************************************************************************
Neural network training  using  modified  Levenberg-Marquardt  with  exact
Hessian calculation and regularization. Subroutine trains  neural  network
with restarts from random positions. Algorithm is well  suited  for  small
and medium scale problems (hundreds of weights).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -9, if internal matrix inverse subroutine failed
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlptrainlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);


/*************************************************************************
Neural  network  training  using  L-BFGS  algorithm  with  regularization.
Subroutine  trains  neural  network  with  restarts from random positions.
Algorithm  is  well  suited  for  problems  of  any dimensionality (memory
requirements and step complexity are linear by weights number).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.
    WStep       -   stopping criterion. Algorithm stops if  step  size  is
                    less than WStep. Recommended value - 0.01.  Zero  step
                    size means stopping after MaxIts iterations.
    MaxIts      -   stopping   criterion.  Algorithm  stops  after  MaxIts
                    iterations (NOT gradient  calculations).  Zero  MaxIts
                    means stopping when step is sufficiently small.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlptrainlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep);


/*************************************************************************
Neural network training using early stopping (base algorithm - L-BFGS with
regularization).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    TrnXY       -   training set
    TrnSize     -   training set size, TrnSize>0
    ValXY       -   validation set
    ValSize     -   validation set size, ValSize>0
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts, either:
                    * strictly positive number - algorithm make specified
                      number of restarts from random position.
                    * -1, in which case algorithm makes exactly one run
                      from the initial state of the network (no randomization).
                    If you don't know what Restarts to choose, choose one
                    one the following:
                    * -1 (deterministic start)
                    * +1 (one random restart)
                    * +5 (moderate amount of random restarts)

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1, ...).
                    *  2, task has been solved, stopping  criterion  met -
                          sufficiently small step size.  Not expected  (we
                          use  EARLY  stopping)  but  possible  and not an
                          error.
                    *  6, task has been solved, stopping  criterion  met -
                          increasing of validation set error.
    Rep         -   training report

NOTE:

Algorithm stops if validation set error increases for  a  long  enough  or
step size is small enought  (there  are  task  where  validation  set  may
decrease for eternity). In any case solution returned corresponds  to  the
minimum of validation set error.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlptraines(const multilayerperceptron &network, const real_2d_array &trnxy, const ae_int_t trnsize, const real_2d_array &valxy, const ae_int_t valsize, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);


/*************************************************************************
Cross-validation estimate of generalization error.

Base algorithm - L-BFGS.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcvlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep);


/*************************************************************************
Cross-validation estimate of generalization error.

Base algorithm - Levenberg-Marquardt.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcvlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep);


/*************************************************************************
This function estimates generalization error using cross-validation on the
current dataset with current training settings.

FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support (C++ computational core)
  !
  ! Second improvement gives constant  speedup (2-3X).  First  improvement
  ! gives  close-to-linear  speedup  on   multicore   systems.   Following
  ! operations can be executed in parallel:
  ! * FoldsCount cross-validation rounds (always)
  ! * NRestarts training sessions performed within each of
  !   cross-validation rounds (if NRestarts>1)
  ! * gradient calculation over large dataset (if dataset is large enough)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.

INPUT PARAMETERS:
    S           -   trainer object
    Network     -   neural network. It must have same number of inputs and
                    output/classes as was specified during creation of the
                    trainer object. Network is not changed  during  cross-
                    validation and is not trained - it  is  used  only  as
                    representative of its architecture. I.e., we  estimate
                    generalization properties of  ARCHITECTURE,  not  some
                    specific network.
    NRestarts   -   number of restarts, >=0:
                    * NRestarts>0  means  that  for  each cross-validation
                      round   specified  number   of  random  restarts  is
                      performed,  with  best  network  being  chosen after
                      training.
                    * NRestarts=0 is same as NRestarts=1
    FoldsCount  -   number of folds in k-fold cross-validation:
                    * 2<=FoldsCount<=size of dataset
                    * recommended value: 10.
                    * values larger than dataset size will be silently
                      truncated down to dataset size

OUTPUT PARAMETERS:
    Rep         -   structure which contains cross-validation estimates:
                    * Rep.RelCLSError - fraction of misclassified cases.
                    * Rep.AvgCE - acerage cross-entropy
                    * Rep.RMSError - root-mean-square error
                    * Rep.AvgError - average error
                    * Rep.AvgRelError - average relative error

NOTE: when no dataset was specified with MLPSetDataset/SetSparseDataset(),
      or subset with only one point  was  given,  zeros  are  returned  as
      estimates.

NOTE: this method performs FoldsCount cross-validation  rounds,  each  one
      with NRestarts random starts.  Thus,  FoldsCount*NRestarts  networks
      are trained in total.

NOTE: Rep.RelCLSError/Rep.AvgCE are zero on regression problems.

NOTE: on classification problems Rep.RMSError/Rep.AvgError/Rep.AvgRelError
      contain errors in prediction of posterior probabilities.

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcv(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, const ae_int_t foldscount, mlpreport &rep);
void smp_mlpkfoldcv(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, const ae_int_t foldscount, mlpreport &rep);


/*************************************************************************
Creation of the network trainer object for regression networks

INPUT PARAMETERS:
    NIn         -   number of inputs, NIn>=1
    NOut        -   number of outputs, NOut>=1

OUTPUT PARAMETERS:
    S           -   neural network trainer object.
                    This structure can be used to train any regression
                    network with NIn inputs and NOut outputs.

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpcreatetrainer(const ae_int_t nin, const ae_int_t nout, mlptrainer &s);


/*************************************************************************
Creation of the network trainer object for classification networks

INPUT PARAMETERS:
    NIn         -   number of inputs, NIn>=1
    NClasses    -   number of classes, NClasses>=2

OUTPUT PARAMETERS:
    S           -   neural network trainer object.
                    This structure can be used to train any classification
                    network with NIn inputs and NOut outputs.

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpcreatetrainercls(const ae_int_t nin, const ae_int_t nclasses, mlptrainer &s);


/*************************************************************************
This function sets "current dataset" of the trainer object to  one  passed
by user.

INPUT PARAMETERS:
    S           -   trainer object
    XY          -   training  set,  see  below  for  information  on   the
                    training set format. This function checks  correctness
                    of  the  dataset  (no  NANs/INFs,  class  numbers  are
                    correct) and throws exception when  incorrect  dataset
                    is passed.
    NPoints     -   points count, >=0.

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
datasetformat is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpsetdataset(const mlptrainer &s, const real_2d_array &xy, const ae_int_t npoints);


/*************************************************************************
This function sets "current dataset" of the trainer object to  one  passed
by user (sparse matrix is used to store dataset).

INPUT PARAMETERS:
    S           -   trainer object
    XY          -   training  set,  see  below  for  information  on   the
                    training set format. This function checks  correctness
                    of  the  dataset  (no  NANs/INFs,  class  numbers  are
                    correct) and throws exception when  incorrect  dataset
                    is passed. Any  sparse  storage  format  can be  used:
                    Hash-table, CRS...
    NPoints     -   points count, >=0

DATASET FORMAT:

This  function  uses  two  different  dataset formats - one for regression
networks, another one for classification networks.

For regression networks with NIn inputs and NOut outputs following dataset
format is used:
* dataset is given by NPoints*(NIn+NOut) matrix
* each row corresponds to one example
* first NIn columns are inputs, next NOut columns are outputs

For classification networks with NIn inputs and NClasses clases  following
datasetformat is used:
* dataset is given by NPoints*(NIn+1) matrix
* each row corresponds to one example
* first NIn columns are inputs, last column stores class number (from 0 to
  NClasses-1).

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpsetsparsedataset(const mlptrainer &s, const sparsematrix &xy, const ae_int_t npoints);


/*************************************************************************
This function sets weight decay coefficient which is used for training.

INPUT PARAMETERS:
    S           -   trainer object
    Decay       -   weight  decay  coefficient,  >=0.  Weight  decay  term
                    'Decay*||Weights||^2' is added to error  function.  If
                    you don't know what Decay to choose, use 1.0E-3.
                    Weight decay can be set to zero,  in this case network
                    is trained without weight decay.

NOTE: by default network uses some small nonzero value for weight decay.

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpsetdecay(const mlptrainer &s, const double decay);


/*************************************************************************
This function sets stopping criteria for the optimizer.

INPUT PARAMETERS:
    S           -   trainer object
    WStep       -   stopping criterion. Algorithm stops if  step  size  is
                    less than WStep. Recommended value - 0.01.  Zero  step
                    size means stopping after MaxIts iterations.
                    WStep>=0.
    MaxIts      -   stopping   criterion.  Algorithm  stops  after  MaxIts
                    epochs (full passes over entire dataset).  Zero MaxIts
                    means stopping when step is sufficiently small.
                    MaxIts>=0.

NOTE: by default, WStep=0.005 and MaxIts=0 are used. These values are also
      used when MLPSetCond() is called with WStep=0 and MaxIts=0.

NOTE: these stopping criteria are used for all kinds of neural training  -
      from "conventional" networks to early stopping ensembles. When  used
      for "conventional" networks, they are  used  as  the  only  stopping
      criteria. When combined with early stopping, they used as ADDITIONAL
      stopping criteria which can terminate early stopping algorithm.

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpsetcond(const mlptrainer &s, const double wstep, const ae_int_t maxits);


/*************************************************************************
This function sets training algorithm: batch training using L-BFGS will be
used.

This algorithm:
* the most robust for small-scale problems, but may be too slow for  large
  scale ones.
* perfoms full pass through the dataset before performing step
* uses conditions specified by MLPSetCond() for stopping
* is default one used by trainer object

INPUT PARAMETERS:
    S           -   trainer object

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpsetalgobatch(const mlptrainer &s);


/*************************************************************************
This function trains neural network passed to this function, using current
dataset (one which was passed to MLPSetDataset() or MLPSetSparseDataset())
and current training settings. Training  from  NRestarts  random  starting
positions is performed, best network is chosen.

Training is performed using current training algorithm.

FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support (C++ computational core)
  !
  ! Second improvement gives constant  speedup (2-3X).  First  improvement
  ! gives  close-to-linear  speedup  on   multicore   systems.   Following
  ! operations can be executed in parallel:
  ! * NRestarts training sessions performed within each of
  !   cross-validation rounds (if NRestarts>1)
  ! * gradient calculation over large dataset (if dataset is large enough)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.

INPUT PARAMETERS:
    S           -   trainer object
    Network     -   neural network. It must have same number of inputs and
                    output/classes as was specified during creation of the
                    trainer object.
    NRestarts   -   number of restarts, >=0:
                    * NRestarts>0 means that specified  number  of  random
                      restarts are performed, best network is chosen after
                      training
                    * NRestarts=0 means that current state of the  network
                      is used for training.

OUTPUT PARAMETERS:
    Network     -   trained network

NOTE: when no dataset was specified with MLPSetDataset/SetSparseDataset(),
      network  is  filled  by zero  values.  Same  behavior  for functions
      MLPStartTraining and MLPContinueTraining.

NOTE: this method uses sum-of-squares error function for training.

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlptrainnetwork(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, mlpreport &rep);
void smp_mlptrainnetwork(const mlptrainer &s, const multilayerperceptron &network, const ae_int_t nrestarts, mlpreport &rep);


/*************************************************************************
IMPORTANT: this is an "expert" version of the MLPTrain() function.  We  do
           not recommend you to use it unless you are pretty sure that you
           need ability to monitor training progress.

This function performs step-by-step training of the neural  network.  Here
"step-by-step" means that training  starts  with  MLPStartTraining() call,
and then user subsequently calls MLPContinueTraining() to perform one more
iteration of the training.

After call to this function trainer object remembers network and  is ready
to  train  it.  However,  no  training  is  performed  until first call to
MLPContinueTraining() function. Subsequent calls  to MLPContinueTraining()
will advance training progress one iteration further.

EXAMPLE:
    >
    > ...initialize network and trainer object....
    >
    > MLPStartTraining(Trainer, Network, True)
    > while MLPContinueTraining(Trainer, Network) do
    >     ...visualize training progress...
    >

INPUT PARAMETERS:
    S           -   trainer object
    Network     -   neural network. It must have same number of inputs and
                    output/classes as was specified during creation of the
                    trainer object.
    RandomStart -   randomize network before training or not:
                    * True  means  that  network  is  randomized  and  its
                      initial state (one which was passed to  the  trainer
                      object) is lost.
                    * False  means  that  training  is  started  from  the
                      current state of the network

OUTPUT PARAMETERS:
    Network     -   neural network which is ready to training (weights are
                    initialized, preprocessor is initialized using current
                    training set)

NOTE: this method uses sum-of-squares error function for training.

NOTE: it is expected that trainer object settings are NOT  changed  during
      step-by-step training, i.e. no  one  changes  stopping  criteria  or
      training set during training. It is possible and there is no defense
      against  such  actions,  but  algorithm  behavior  in  such cases is
      undefined and can be unpredictable.

  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
void mlpstarttraining(const mlptrainer &s, const multilayerperceptron &network, const bool randomstart);


/*************************************************************************
IMPORTANT: this is an "expert" version of the MLPTrain() function.  We  do
           not recommend you to use it unless you are pretty sure that you
           need ability to monitor training progress.

FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support (C++ computational core)
  !
  ! Second improvement gives constant  speedup (2-3X).  First  improvement
  ! gives  close-to-linear  speedup  on   multicore   systems.   Following
  ! operations can be executed in parallel:
  ! * gradient calculation over large dataset (if dataset is large enough)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.

This function performs step-by-step training of the neural  network.  Here
"step-by-step" means that training starts  with  MLPStartTraining()  call,
and then user subsequently calls MLPContinueTraining() to perform one more
iteration of the training.

This  function  performs  one  more  iteration of the training and returns
either True (training continues) or False (training stopped). In case True
was returned, Network weights are updated according to the  current  state
of the optimization progress. In case False was  returned,  no  additional
updates is performed (previous update of  the  network weights moved us to
the final point, and no additional updates is needed).

EXAMPLE:
    >
    > [initialize network and trainer object]
    >
    > MLPStartTraining(Trainer, Network, True)
    > while MLPContinueTraining(Trainer, Network) do
    >     [visualize training progress]
    >

INPUT PARAMETERS:
    S           -   trainer object
    Network     -   neural  network  structure,  which  is  used to  store
                    current state of the training process.

OUTPUT PARAMETERS:
    Network     -   weights of the neural network  are  rewritten  by  the
                    current approximation.

NOTE: this method uses sum-of-squares error function for training.

NOTE: it is expected that trainer object settings are NOT  changed  during
      step-by-step training, i.e. no  one  changes  stopping  criteria  or
      training set during training. It is possible and there is no defense
      against  such  actions,  but  algorithm  behavior  in  such cases is
      undefined and can be unpredictable.

NOTE: It  is  expected that Network is the same one which  was  passed  to
      MLPStartTraining() function.  However,  THIS  function  checks  only
      following:
      * that number of network inputs is consistent with trainer object
        settings
      * that number of network outputs/classes is consistent with  trainer
        object settings
      * that number of network weights is the same as number of weights in
        the network passed to MLPStartTraining() function
      Exception is thrown when these conditions are violated.

      It is also expected that you do not change state of the  network  on
      your own - the only party who has right to change network during its
      training is a trainer object. Any attempt to interfere with  trainer
      may lead to unpredictable results.


  -- ALGLIB --
     Copyright 23.07.2012 by Bochkanov Sergey
*************************************************************************/
bool mlpcontinuetraining(const mlptrainer &s, const multilayerperceptron &network);
bool smp_mlpcontinuetraining(const mlptrainer &s, const multilayerperceptron &network);


/*************************************************************************
Training neural networks ensemble using  bootstrap  aggregating (bagging).
Modified Levenberg-Marquardt algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpebagginglm(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors);


/*************************************************************************
Training neural networks ensemble using  bootstrap  aggregating (bagging).
L-BFGS algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpebagginglbfgs(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors);


/*************************************************************************
Training neural networks ensemble using early stopping.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  6, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlpetraines(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);


/*************************************************************************
This function trains neural network ensemble passed to this function using
current dataset and early stopping training algorithm. Each early stopping
round performs NRestarts  random  restarts  (thus,  EnsembleSize*NRestarts
training rounds is performed in total).

FOR USERS OF COMMERCIAL EDITION:

  ! Commercial version of ALGLIB includes two  important  improvements  of
  ! this function:
  ! * multicore support (C++ and C# computational cores)
  ! * SSE support (C++ computational core)
  !
  ! Second improvement gives constant  speedup (2-3X).  First  improvement
  ! gives  close-to-linear  speedup  on   multicore   systems.   Following
  ! operations can be executed in parallel:
  ! * EnsembleSize  training  sessions  performed  for  each  of  ensemble
  !   members (always parallelized)
  ! * NRestarts  training  sessions  performed  within  each  of  training
  !   sessions (if NRestarts>1)
  ! * gradient calculation over large dataset (if dataset is large enough)
  !
  ! In order to use multicore features you have to:
  ! * use commercial version of ALGLIB
  ! * call  this  function  with  "smp_"  prefix,  which  indicates  that
  !   multicore code will be used (for multicore support)
  !
  ! In order to use SSE features you have to:
  ! * use commercial version of ALGLIB on Intel processors
  ! * use C++ computational core
  !
  ! This note is given for users of commercial edition; if  you  use  GPL
  ! edition, you still will be able to call smp-version of this function,
  ! but all computations will be done serially.
  !
  ! We recommend you to carefully read ALGLIB Reference  Manual,  section
  ! called 'SMP support', before using parallel version of this function.

INPUT PARAMETERS:
    S           -   trainer object;
    Ensemble    -   neural network ensemble. It must have same  number  of
                    inputs and outputs/classes  as  was  specified  during
                    creation of the trainer object.
    NRestarts   -   number of restarts, >=0:
                    * NRestarts>0 means that specified  number  of  random
                      restarts are performed during each ES round;
                    * NRestarts=0 is silently replaced by 1.

OUTPUT PARAMETERS:
    Ensemble    -   trained ensemble;
    Rep         -   it contains all type of errors.

NOTE: this training method uses BOTH early stopping and weight decay!  So,
      you should select weight decay before starting training just as  you
      select it before training "conventional" networks.

NOTE: when no dataset was specified with MLPSetDataset/SetSparseDataset(),
      or  single-point  dataset  was  passed,  ensemble  is filled by zero
      values.

NOTE: this method uses sum-of-squares error function for training.

  -- ALGLIB --
     Copyright 22.08.2012 by Bochkanov Sergey
*************************************************************************/
void mlptrainensemblees(const mlptrainer &s, const mlpensemble &ensemble, const ae_int_t nrestarts, mlpreport &rep);
void smp_mlptrainensemblees(const mlptrainer &s, const mlpensemble &ensemble, const ae_int_t nrestarts, mlpreport &rep);

/*************************************************************************
Principal components analysis

Subroutine  builds  orthogonal  basis  where  first  axis  corresponds  to
direction with maximum variance, second axis maximizes variance in subspace
orthogonal to first axis and so on.

It should be noted that, unlike LDA, PCA does not use class labels.

COMMERCIAL EDITION OF ALGLIB:

  ! Commercial version of ALGLIB includes one  important  improvement   of
  ! this function, which can be used from C++ and C#:
  ! * Intel MKL support (lightweight Intel MKL is shipped with ALGLIB)
  !
  ! Intel MKL gives approximately constant  (with  respect  to  number  of
  ! worker threads) acceleration factor which depends on CPU  being  used,
  ! problem  size  and  "baseline"  ALGLIB  edition  which  is  used   for
  ! comparison. Best results are achieved  for  high-dimensional  problems
  ! (NVars is at least 256).
  !
  ! We recommend you to read 'Working with commercial version' section  of
  ! ALGLIB Reference Manual in order to find out how to  use  performance-
  ! related features provided by commercial edition of ALGLIB.

INPUT PARAMETERS:
    X           -   dataset, array[0..NPoints-1,0..NVars-1].
                    matrix contains ONLY INDEPENDENT VARIABLES.
    NPoints     -   dataset size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if SVD subroutine haven't converged
                    * -1, if wrong parameters has been passed (NPoints<0,
                          NVars<1)
                    *  1, if task is solved
    S2          -   array[0..NVars-1]. variance values corresponding
                    to basis vectors.
    V           -   array[0..NVars-1,0..NVars-1]
                    matrix, whose columns store basis vectors.

  -- ALGLIB --
     Copyright 25.08.2008 by Bochkanov Sergey
*************************************************************************/
void pcabuildbasis(const real_2d_array &x, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, real_1d_array &s2, real_2d_array &v);
}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (FUNCTIONS)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
void dserrallocate(ae_int_t nclasses,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
void dserraccumulate(/* Real    */ ae_vector* buf,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* desiredy,
     ae_state *_state);
void dserrfinish(/* Real    */ ae_vector* buf, ae_state *_state);
void dsnormalize(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     /* Real    */ ae_vector* means,
     /* Real    */ ae_vector* sigmas,
     ae_state *_state);
void dsnormalizec(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     /* Real    */ ae_vector* means,
     /* Real    */ ae_vector* sigmas,
     ae_state *_state);
double dsgetmeanmindistance(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_state *_state);
void dstie(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Integer */ ae_vector* ties,
     ae_int_t* tiecount,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     ae_state *_state);
void dstiefasti(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t n,
     /* Integer */ ae_vector* ties,
     ae_int_t* tiecount,
     /* Real    */ ae_vector* bufr,
     /* Integer */ ae_vector* bufi,
     ae_state *_state);
void dsoptimalsplit2(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* c,
     ae_int_t n,
     ae_int_t* info,
     double* threshold,
     double* pal,
     double* pbl,
     double* par,
     double* pbr,
     double* cve,
     ae_state *_state);
void dsoptimalsplit2fast(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* c,
     /* Integer */ ae_vector* tiesbuf,
     /* Integer */ ae_vector* cntbuf,
     /* Real    */ ae_vector* bufr,
     /* Integer */ ae_vector* bufi,
     ae_int_t n,
     ae_int_t nc,
     double alpha,
     ae_int_t* info,
     double* threshold,
     double* rms,
     double* cvrms,
     ae_state *_state);
void dssplitk(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* c,
     ae_int_t n,
     ae_int_t nc,
     ae_int_t kmax,
     ae_int_t* info,
     /* Real    */ ae_vector* thresholds,
     ae_int_t* ni,
     double* cve,
     ae_state *_state);
void dsoptimalsplitk(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* c,
     ae_int_t n,
     ae_int_t nc,
     ae_int_t kmax,
     ae_int_t* info,
     /* Real    */ ae_vector* thresholds,
     ae_int_t* ni,
     double* cve,
     ae_state *_state);
void _cvreport_init(void* _p, ae_state *_state);
void _cvreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _cvreport_clear(void* _p);
void _cvreport_destroy(void* _p);
void clusterizercreate(clusterizerstate* s, ae_state *_state);
void clusterizersetpoints(clusterizerstate* s,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nfeatures,
     ae_int_t disttype,
     ae_state *_state);
void clusterizersetdistances(clusterizerstate* s,
     /* Real    */ ae_matrix* d,
     ae_int_t npoints,
     ae_bool isupper,
     ae_state *_state);
void clusterizersetahcalgo(clusterizerstate* s,
     ae_int_t algo,
     ae_state *_state);
void clusterizersetkmeanslimits(clusterizerstate* s,
     ae_int_t restarts,
     ae_int_t maxits,
     ae_state *_state);
void clusterizersetkmeansinit(clusterizerstate* s,
     ae_int_t initalgo,
     ae_state *_state);
void clusterizerrunahc(clusterizerstate* s,
     ahcreport* rep,
     ae_state *_state);
void _pexec_clusterizerrunahc(clusterizerstate* s,
    ahcreport* rep, ae_state *_state);
void clusterizerrunkmeans(clusterizerstate* s,
     ae_int_t k,
     kmeansreport* rep,
     ae_state *_state);
void _pexec_clusterizerrunkmeans(clusterizerstate* s,
    ae_int_t k,
    kmeansreport* rep, ae_state *_state);
void clusterizergetdistances(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nfeatures,
     ae_int_t disttype,
     /* Real    */ ae_matrix* d,
     ae_state *_state);
void _pexec_clusterizergetdistances(/* Real    */ ae_matrix* xy,
    ae_int_t npoints,
    ae_int_t nfeatures,
    ae_int_t disttype,
    /* Real    */ ae_matrix* d, ae_state *_state);
void clusterizergetdistancesbuf(apbuffers* buf,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nfeatures,
     ae_int_t disttype,
     /* Real    */ ae_matrix* d,
     ae_state *_state);
void clusterizergetkclusters(ahcreport* rep,
     ae_int_t k,
     /* Integer */ ae_vector* cidx,
     /* Integer */ ae_vector* cz,
     ae_state *_state);
void clusterizerseparatedbydist(ahcreport* rep,
     double r,
     ae_int_t* k,
     /* Integer */ ae_vector* cidx,
     /* Integer */ ae_vector* cz,
     ae_state *_state);
void clusterizerseparatedbycorr(ahcreport* rep,
     double r,
     ae_int_t* k,
     /* Integer */ ae_vector* cidx,
     /* Integer */ ae_vector* cz,
     ae_state *_state);
void kmeansinitbuf(kmeansbuffers* buf, ae_state *_state);
void kmeansgenerateinternal(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t k,
     ae_int_t initalgo,
     ae_int_t maxits,
     ae_int_t restarts,
     ae_bool kmeansdbgnoits,
     ae_int_t* info,
     ae_int_t* iterationscount,
     /* Real    */ ae_matrix* ccol,
     ae_bool needccol,
     /* Real    */ ae_matrix* crow,
     ae_bool needcrow,
     /* Integer */ ae_vector* xyc,
     double* energy,
     kmeansbuffers* buf,
     ae_state *_state);
void kmeansupdatedistances(/* Real    */ ae_matrix* xy,
     ae_int_t idx0,
     ae_int_t idx1,
     ae_int_t nvars,
     /* Real    */ ae_matrix* ct,
     ae_int_t cidx0,
     ae_int_t cidx1,
     /* Integer */ ae_vector* xyc,
     /* Real    */ ae_vector* xydist2,
     ae_shared_pool* bufferpool,
     ae_state *_state);
void _kmeansbuffers_init(void* _p, ae_state *_state);
void _kmeansbuffers_init_copy(void* _dst, void* _src, ae_state *_state);
void _kmeansbuffers_clear(void* _p);
void _kmeansbuffers_destroy(void* _p);
void _clusterizerstate_init(void* _p, ae_state *_state);
void _clusterizerstate_init_copy(void* _dst, void* _src, ae_state *_state);
void _clusterizerstate_clear(void* _p);
void _clusterizerstate_destroy(void* _p);
void _ahcreport_init(void* _p, ae_state *_state);
void _ahcreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _ahcreport_clear(void* _p);
void _ahcreport_destroy(void* _p);
void _kmeansreport_init(void* _p, ae_state *_state);
void _kmeansreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _kmeansreport_clear(void* _p);
void _kmeansreport_destroy(void* _p);
void kmeansgenerate(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t k,
     ae_int_t restarts,
     ae_int_t* info,
     /* Real    */ ae_matrix* c,
     /* Integer */ ae_vector* xyc,
     ae_state *_state);
void dfbuildrandomdecisionforest(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t ntrees,
     double r,
     ae_int_t* info,
     decisionforest* df,
     dfreport* rep,
     ae_state *_state);
void dfbuildrandomdecisionforestx1(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t ntrees,
     ae_int_t nrndvars,
     double r,
     ae_int_t* info,
     decisionforest* df,
     dfreport* rep,
     ae_state *_state);
void dfbuildinternal(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t ntrees,
     ae_int_t samplesize,
     ae_int_t nfeatures,
     ae_int_t flags,
     ae_int_t* info,
     decisionforest* df,
     dfreport* rep,
     ae_state *_state);
void dfprocess(decisionforest* df,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void dfprocessi(decisionforest* df,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
double dfrelclserror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double dfavgce(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double dfrmserror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double dfavgerror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double dfavgrelerror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
void dfcopy(decisionforest* df1, decisionforest* df2, ae_state *_state);
void dfalloc(ae_serializer* s, decisionforest* forest, ae_state *_state);
void dfserialize(ae_serializer* s,
     decisionforest* forest,
     ae_state *_state);
void dfunserialize(ae_serializer* s,
     decisionforest* forest,
     ae_state *_state);
void _decisionforest_init(void* _p, ae_state *_state);
void _decisionforest_init_copy(void* _dst, void* _src, ae_state *_state);
void _decisionforest_clear(void* _p);
void _decisionforest_destroy(void* _p);
void _dfreport_init(void* _p, ae_state *_state);
void _dfreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _dfreport_clear(void* _p);
void _dfreport_destroy(void* _p);
void _dfinternalbuffers_init(void* _p, ae_state *_state);
void _dfinternalbuffers_init_copy(void* _dst, void* _src, ae_state *_state);
void _dfinternalbuffers_clear(void* _p);
void _dfinternalbuffers_destroy(void* _p);
void lrbuild(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state);
void lrbuilds(/* Real    */ ae_matrix* xy,
     /* Real    */ ae_vector* s,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state);
void lrbuildzs(/* Real    */ ae_matrix* xy,
     /* Real    */ ae_vector* s,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state);
void lrbuildz(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state);
void lrunpack(linearmodel* lm,
     /* Real    */ ae_vector* v,
     ae_int_t* nvars,
     ae_state *_state);
void lrpack(/* Real    */ ae_vector* v,
     ae_int_t nvars,
     linearmodel* lm,
     ae_state *_state);
double lrprocess(linearmodel* lm,
     /* Real    */ ae_vector* x,
     ae_state *_state);
double lrrmserror(linearmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double lravgerror(linearmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double lravgrelerror(linearmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
void lrcopy(linearmodel* lm1, linearmodel* lm2, ae_state *_state);
void lrlines(/* Real    */ ae_matrix* xy,
     /* Real    */ ae_vector* s,
     ae_int_t n,
     ae_int_t* info,
     double* a,
     double* b,
     double* vara,
     double* varb,
     double* covab,
     double* corrab,
     double* p,
     ae_state *_state);
void lrline(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t* info,
     double* a,
     double* b,
     ae_state *_state);
void _linearmodel_init(void* _p, ae_state *_state);
void _linearmodel_init_copy(void* _dst, void* _src, ae_state *_state);
void _linearmodel_clear(void* _p);
void _linearmodel_destroy(void* _p);
void _lrreport_init(void* _p, ae_state *_state);
void _lrreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _lrreport_clear(void* _p);
void _lrreport_destroy(void* _p);
void filtersma(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_int_t k,
     ae_state *_state);
void filterema(/* Real    */ ae_vector* x,
     ae_int_t n,
     double alpha,
     ae_state *_state);
void filterlrma(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_int_t k,
     ae_state *_state);
void fisherlda(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     /* Real    */ ae_vector* w,
     ae_state *_state);
void fisherldan(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     /* Real    */ ae_matrix* w,
     ae_state *_state);
void _pexec_fisherldan(/* Real    */ ae_matrix* xy,
    ae_int_t npoints,
    ae_int_t nvars,
    ae_int_t nclasses,
    ae_int_t* info,
    /* Real    */ ae_matrix* w, ae_state *_state);
ae_int_t mlpgradsplitcost(ae_state *_state);
ae_int_t mlpgradsplitsize(ae_state *_state);
void mlpcreate0(ae_int_t nin,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreate1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreate2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreateb0(ae_int_t nin,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreateb1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreateb2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreater0(ae_int_t nin,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreater1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreater2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreatec0(ae_int_t nin,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreatec1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreatec2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcopy(multilayerperceptron* network1,
     multilayerperceptron* network2,
     ae_state *_state);
void mlpcopyshared(multilayerperceptron* network1,
     multilayerperceptron* network2,
     ae_state *_state);
ae_bool mlpsamearchitecture(multilayerperceptron* network1,
     multilayerperceptron* network2,
     ae_state *_state);
void mlpcopytunableparameters(multilayerperceptron* network1,
     multilayerperceptron* network2,
     ae_state *_state);
void mlpexporttunableparameters(multilayerperceptron* network,
     /* Real    */ ae_vector* p,
     ae_int_t* pcount,
     ae_state *_state);
void mlpimporttunableparameters(multilayerperceptron* network,
     /* Real    */ ae_vector* p,
     ae_state *_state);
void mlpserializeold(multilayerperceptron* network,
     /* Real    */ ae_vector* ra,
     ae_int_t* rlen,
     ae_state *_state);
void mlpunserializeold(/* Real    */ ae_vector* ra,
     multilayerperceptron* network,
     ae_state *_state);
void mlprandomize(multilayerperceptron* network, ae_state *_state);
void mlprandomizefull(multilayerperceptron* network, ae_state *_state);
void mlpinitpreprocessor(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state);
void mlpinitpreprocessorsparse(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t ssize,
     ae_state *_state);
void mlpinitpreprocessorsubset(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t setsize,
     /* Integer */ ae_vector* idx,
     ae_int_t subsetsize,
     ae_state *_state);
void mlpinitpreprocessorsparsesubset(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t setsize,
     /* Integer */ ae_vector* idx,
     ae_int_t subsetsize,
     ae_state *_state);
void mlpproperties(multilayerperceptron* network,
     ae_int_t* nin,
     ae_int_t* nout,
     ae_int_t* wcount,
     ae_state *_state);
ae_int_t mlpntotal(multilayerperceptron* network, ae_state *_state);
ae_int_t mlpgetinputscount(multilayerperceptron* network,
     ae_state *_state);
ae_int_t mlpgetoutputscount(multilayerperceptron* network,
     ae_state *_state);
ae_int_t mlpgetweightscount(multilayerperceptron* network,
     ae_state *_state);
ae_bool mlpissoftmax(multilayerperceptron* network, ae_state *_state);
ae_int_t mlpgetlayerscount(multilayerperceptron* network,
     ae_state *_state);
ae_int_t mlpgetlayersize(multilayerperceptron* network,
     ae_int_t k,
     ae_state *_state);
void mlpgetinputscaling(multilayerperceptron* network,
     ae_int_t i,
     double* mean,
     double* sigma,
     ae_state *_state);
void mlpgetoutputscaling(multilayerperceptron* network,
     ae_int_t i,
     double* mean,
     double* sigma,
     ae_state *_state);
void mlpgetneuroninfo(multilayerperceptron* network,
     ae_int_t k,
     ae_int_t i,
     ae_int_t* fkind,
     double* threshold,
     ae_state *_state);
double mlpgetweight(multilayerperceptron* network,
     ae_int_t k0,
     ae_int_t i0,
     ae_int_t k1,
     ae_int_t i1,
     ae_state *_state);
void mlpsetinputscaling(multilayerperceptron* network,
     ae_int_t i,
     double mean,
     double sigma,
     ae_state *_state);
void mlpsetoutputscaling(multilayerperceptron* network,
     ae_int_t i,
     double mean,
     double sigma,
     ae_state *_state);
void mlpsetneuroninfo(multilayerperceptron* network,
     ae_int_t k,
     ae_int_t i,
     ae_int_t fkind,
     double threshold,
     ae_state *_state);
void mlpsetweight(multilayerperceptron* network,
     ae_int_t k0,
     ae_int_t i0,
     ae_int_t k1,
     ae_int_t i1,
     double w,
     ae_state *_state);
void mlpactivationfunction(double net,
     ae_int_t k,
     double* f,
     double* df,
     double* d2f,
     ae_state *_state);
void mlpprocess(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mlpprocessi(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
double mlperror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlperror(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlperrorsparse(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlperrorsparse(multilayerperceptron* network,
    sparsematrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlperrorn(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state);
ae_int_t mlpclserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
ae_int_t _pexec_mlpclserror(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlprelclserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlprelclserror(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlprelclserrorsparse(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlprelclserrorsparse(multilayerperceptron* network,
    sparsematrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlpavgce(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlpavgce(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlpavgcesparse(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlpavgcesparse(multilayerperceptron* network,
    sparsematrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlprmserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlprmserror(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlprmserrorsparse(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlprmserrorsparse(multilayerperceptron* network,
    sparsematrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlpavgerror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlpavgerror(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlpavgerrorsparse(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlpavgerrorsparse(multilayerperceptron* network,
    sparsematrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlpavgrelerror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlpavgrelerror(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t npoints, ae_state *_state);
double mlpavgrelerrorsparse(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double _pexec_mlpavgrelerrorsparse(multilayerperceptron* network,
    sparsematrix* xy,
    ae_int_t npoints, ae_state *_state);
void mlpgrad(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* desiredy,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void mlpgradn(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* desiredy,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void mlpgradbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void _pexec_mlpgradbatch(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t ssize,
    double* e,
    /* Real    */ ae_vector* grad, ae_state *_state);
void mlpgradbatchsparse(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void _pexec_mlpgradbatchsparse(multilayerperceptron* network,
    sparsematrix* xy,
    ae_int_t ssize,
    double* e,
    /* Real    */ ae_vector* grad, ae_state *_state);
void mlpgradbatchsubset(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t setsize,
     /* Integer */ ae_vector* idx,
     ae_int_t subsetsize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void _pexec_mlpgradbatchsubset(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t setsize,
    /* Integer */ ae_vector* idx,
    ae_int_t subsetsize,
    double* e,
    /* Real    */ ae_vector* grad, ae_state *_state);
void mlpgradbatchsparsesubset(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t setsize,
     /* Integer */ ae_vector* idx,
     ae_int_t subsetsize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void _pexec_mlpgradbatchsparsesubset(multilayerperceptron* network,
    sparsematrix* xy,
    ae_int_t setsize,
    /* Integer */ ae_vector* idx,
    ae_int_t subsetsize,
    double* e,
    /* Real    */ ae_vector* grad, ae_state *_state);
void mlpgradbatchx(multilayerperceptron* network,
     /* Real    */ ae_matrix* densexy,
     sparsematrix* sparsexy,
     ae_int_t datasetsize,
     ae_int_t datasettype,
     /* Integer */ ae_vector* idx,
     ae_int_t subset0,
     ae_int_t subset1,
     ae_int_t subsettype,
     ae_shared_pool* buf,
     ae_shared_pool* gradbuf,
     ae_state *_state);
void mlpgradnbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void mlphessiannbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state);
void mlphessianbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state);
void mlpinternalprocessvector(/* Integer */ ae_vector* structinfo,
     /* Real    */ ae_vector* weights,
     /* Real    */ ae_vector* columnmeans,
     /* Real    */ ae_vector* columnsigmas,
     /* Real    */ ae_vector* neurons,
     /* Real    */ ae_vector* dfdnet,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mlpalloc(ae_serializer* s,
     multilayerperceptron* network,
     ae_state *_state);
void mlpserialize(ae_serializer* s,
     multilayerperceptron* network,
     ae_state *_state);
void mlpunserialize(ae_serializer* s,
     multilayerperceptron* network,
     ae_state *_state);
void mlpallerrorssubset(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t setsize,
     /* Integer */ ae_vector* subset,
     ae_int_t subsetsize,
     modelerrors* rep,
     ae_state *_state);
void _pexec_mlpallerrorssubset(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t setsize,
    /* Integer */ ae_vector* subset,
    ae_int_t subsetsize,
    modelerrors* rep, ae_state *_state);
void mlpallerrorssparsesubset(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t setsize,
     /* Integer */ ae_vector* subset,
     ae_int_t subsetsize,
     modelerrors* rep,
     ae_state *_state);
void _pexec_mlpallerrorssparsesubset(multilayerperceptron* network,
    sparsematrix* xy,
    ae_int_t setsize,
    /* Integer */ ae_vector* subset,
    ae_int_t subsetsize,
    modelerrors* rep, ae_state *_state);
double mlperrorsubset(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t setsize,
     /* Integer */ ae_vector* subset,
     ae_int_t subsetsize,
     ae_state *_state);
double _pexec_mlperrorsubset(multilayerperceptron* network,
    /* Real    */ ae_matrix* xy,
    ae_int_t setsize,
    /* Integer */ ae_vector* subset,
    ae_int_t subsetsize, ae_state *_state);
double mlperrorsparsesubset(multilayerperceptron* network,
     sparsematrix* xy,
     ae_int_t setsize,
     /* Integer */ ae_vector* subset,
     ae_int_t subsetsize,
     ae_state *_state);
double _pexec_mlperrorsparsesubset(multilayerperceptron* network,
    sparsematrix* xy,
    ae_int_t setsize,
    /* Integer */ ae_vector* subset,
    ae_int_t subsetsize, ae_state *_state);
void mlpallerrorsx(multilayerperceptron* network,
     /* Real    */ ae_matrix* densexy,
     sparsematrix* sparsexy,
     ae_int_t datasetsize,
     ae_int_t datasettype,
     /* Integer */ ae_vector* idx,
     ae_int_t subset0,
     ae_int_t subset1,
     ae_int_t subsettype,
     ae_shared_pool* buf,
     modelerrors* rep,
     ae_state *_state);
void _modelerrors_init(void* _p, ae_state *_state);
void _modelerrors_init_copy(void* _dst, void* _src, ae_state *_state);
void _modelerrors_clear(void* _p);
void _modelerrors_destroy(void* _p);
void _smlpgrad_init(void* _p, ae_state *_state);
void _smlpgrad_init_copy(void* _dst, void* _src, ae_state *_state);
void _smlpgrad_clear(void* _p);
void _smlpgrad_destroy(void* _p);
void _multilayerperceptron_init(void* _p, ae_state *_state);
void _multilayerperceptron_init_copy(void* _dst, void* _src, ae_state *_state);
void _multilayerperceptron_clear(void* _p);
void _multilayerperceptron_destroy(void* _p);
void mnltrainh(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     logitmodel* lm,
     mnlreport* rep,
     ae_state *_state);
void mnlprocess(logitmodel* lm,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mnlprocessi(logitmodel* lm,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mnlunpack(logitmodel* lm,
     /* Real    */ ae_matrix* a,
     ae_int_t* nvars,
     ae_int_t* nclasses,
     ae_state *_state);
void mnlpack(/* Real    */ ae_matrix* a,
     ae_int_t nvars,
     ae_int_t nclasses,
     logitmodel* lm,
     ae_state *_state);
void mnlcopy(logitmodel* lm1, logitmodel* lm2, ae_state *_state);
double mnlavgce(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mnlrelclserror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mnlrmserror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mnlavgerror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mnlavgrelerror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state);
ae_int_t mnlclserror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
void _logitmodel_init(void* _p, ae_state *_state);
void _logitmodel_init_copy(void* _dst, void* _src, ae_state *_state);
void _logitmodel_clear(void* _p);
void _logitmodel_destroy(void* _p);
void _logitmcstate_init(void* _p, ae_state *_state);
void _logitmcstate_init_copy(void* _dst, void* _src, ae_state *_state);
void _logitmcstate_clear(void* _p);
void _logitmcstate_destroy(void* _p);
void _mnlreport_init(void* _p, ae_state *_state);
void _mnlreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _mnlreport_clear(void* _p);
void _mnlreport_destroy(void* _p);
void mcpdcreate(ae_int_t n, mcpdstate* s, ae_state *_state);
void mcpdcreateentry(ae_int_t n,
     ae_int_t entrystate,
     mcpdstate* s,
     ae_state *_state);
void mcpdcreateexit(ae_int_t n,
     ae_int_t exitstate,
     mcpdstate* s,
     ae_state *_state);
void mcpdcreateentryexit(ae_int_t n,
     ae_int_t entrystate,
     ae_int_t exitstate,
     mcpdstate* s,
     ae_state *_state);
void mcpdaddtrack(mcpdstate* s,
     /* Real    */ ae_matrix* xy,
     ae_int_t k,
     ae_state *_state);
void mcpdsetec(mcpdstate* s,
     /* Real    */ ae_matrix* ec,
     ae_state *_state);
void mcpdaddec(mcpdstate* s,
     ae_int_t i,
     ae_int_t j,
     double c,
     ae_state *_state);
void mcpdsetbc(mcpdstate* s,
     /* Real    */ ae_matrix* bndl,
     /* Real    */ ae_matrix* bndu,
     ae_state *_state);
void mcpdaddbc(mcpdstate* s,
     ae_int_t i,
     ae_int_t j,
     double bndl,
     double bndu,
     ae_state *_state);
void mcpdsetlc(mcpdstate* s,
     /* Real    */ ae_matrix* c,
     /* Integer */ ae_vector* ct,
     ae_int_t k,
     ae_state *_state);
void mcpdsettikhonovregularizer(mcpdstate* s, double v, ae_state *_state);
void mcpdsetprior(mcpdstate* s,
     /* Real    */ ae_matrix* pp,
     ae_state *_state);
void mcpdsetpredictionweights(mcpdstate* s,
     /* Real    */ ae_vector* pw,
     ae_state *_state);
void mcpdsolve(mcpdstate* s, ae_state *_state);
void mcpdresults(mcpdstate* s,
     /* Real    */ ae_matrix* p,
     mcpdreport* rep,
     ae_state *_state);
void _mcpdstate_init(void* _p, ae_state *_state);
void _mcpdstate_init_copy(void* _dst, void* _src, ae_state *_state);
void _mcpdstate_clear(void* _p);
void _mcpdstate_destroy(void* _p);
void _mcpdreport_init(void* _p, ae_state *_state);
void _mcpdreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _mcpdreport_clear(void* _p);
void _mcpdreport_destroy(void* _p);
void mlpecreate0(ae_int_t nin,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreate1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreate2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreateb0(ae_int_t nin,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreateb1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreateb2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreater0(ae_int_t nin,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreater1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreater2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreatec0(ae_int_t nin,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreatec1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreatec2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreatefromnetwork(multilayerperceptron* network,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecopy(mlpensemble* ensemble1,
     mlpensemble* ensemble2,
     ae_state *_state);
void mlperandomize(mlpensemble* ensemble, ae_state *_state);
void mlpeproperties(mlpensemble* ensemble,
     ae_int_t* nin,
     ae_int_t* nout,
     ae_state *_state);
ae_bool mlpeissoftmax(mlpensemble* ensemble, ae_state *_state);
void mlpeprocess(mlpensemble* ensemble,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mlpeprocessi(mlpensemble* ensemble,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mlpeallerrorsx(mlpensemble* ensemble,
     /* Real    */ ae_matrix* densexy,
     sparsematrix* sparsexy,
     ae_int_t datasetsize,
     ae_int_t datasettype,
     /* Integer */ ae_vector* idx,
     ae_int_t subset0,
     ae_int_t subset1,
     ae_int_t subsettype,
     ae_shared_pool* buf,
     modelerrors* rep,
     ae_state *_state);
void mlpeallerrorssparse(mlpensemble* ensemble,
     sparsematrix* xy,
     ae_int_t npoints,
     double* relcls,
     double* avgce,
     double* rms,
     double* avg,
     double* avgrel,
     ae_state *_state);
double mlperelclserror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpeavgce(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpermserror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpeavgerror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpeavgrelerror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
void mlpealloc(ae_serializer* s, mlpensemble* ensemble, ae_state *_state);
void mlpeserialize(ae_serializer* s,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpeunserialize(ae_serializer* s,
     mlpensemble* ensemble,
     ae_state *_state);
void _mlpensemble_init(void* _p, ae_state *_state);
void _mlpensemble_init_copy(void* _dst, void* _src, ae_state *_state);
void _mlpensemble_clear(void* _p);
void _mlpensemble_destroy(void* _p);
void mlptrainlm(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state);
void mlptrainlbfgs(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state);
void mlptraines(multilayerperceptron* network,
     /* Real    */ ae_matrix* trnxy,
     ae_int_t trnsize,
     /* Real    */ ae_matrix* valxy,
     ae_int_t valsize,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state);
void mlpkfoldcvlbfgs(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_int_t foldscount,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state);
void mlpkfoldcvlm(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t foldscount,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state);
void mlpkfoldcv(mlptrainer* s,
     multilayerperceptron* network,
     ae_int_t nrestarts,
     ae_int_t foldscount,
     mlpreport* rep,
     ae_state *_state);
void _pexec_mlpkfoldcv(mlptrainer* s,
    multilayerperceptron* network,
    ae_int_t nrestarts,
    ae_int_t foldscount,
    mlpreport* rep, ae_state *_state);
void mlpcreatetrainer(ae_int_t nin,
     ae_int_t nout,
     mlptrainer* s,
     ae_state *_state);
void mlpcreatetrainercls(ae_int_t nin,
     ae_int_t nclasses,
     mlptrainer* s,
     ae_state *_state);
void mlpsetdataset(mlptrainer* s,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
void mlpsetsparsedataset(mlptrainer* s,
     sparsematrix* xy,
     ae_int_t npoints,
     ae_state *_state);
void mlpsetdecay(mlptrainer* s, double decay, ae_state *_state);
void mlpsetcond(mlptrainer* s,
     double wstep,
     ae_int_t maxits,
     ae_state *_state);
void mlpsetalgobatch(mlptrainer* s, ae_state *_state);
void mlptrainnetwork(mlptrainer* s,
     multilayerperceptron* network,
     ae_int_t nrestarts,
     mlpreport* rep,
     ae_state *_state);
void _pexec_mlptrainnetwork(mlptrainer* s,
    multilayerperceptron* network,
    ae_int_t nrestarts,
    mlpreport* rep, ae_state *_state);
void mlpstarttraining(mlptrainer* s,
     multilayerperceptron* network,
     ae_bool randomstart,
     ae_state *_state);
ae_bool mlpcontinuetraining(mlptrainer* s,
     multilayerperceptron* network,
     ae_state *_state);
ae_bool _pexec_mlpcontinuetraining(mlptrainer* s,
    multilayerperceptron* network, ae_state *_state);
void mlpebagginglm(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state);
void mlpebagginglbfgs(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state);
void mlpetraines(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state);
void mlptrainensemblees(mlptrainer* s,
     mlpensemble* ensemble,
     ae_int_t nrestarts,
     mlpreport* rep,
     ae_state *_state);
void _pexec_mlptrainensemblees(mlptrainer* s,
    mlpensemble* ensemble,
    ae_int_t nrestarts,
    mlpreport* rep, ae_state *_state);
void _mlpreport_init(void* _p, ae_state *_state);
void _mlpreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _mlpreport_clear(void* _p);
void _mlpreport_destroy(void* _p);
void _mlpcvreport_init(void* _p, ae_state *_state);
void _mlpcvreport_init_copy(void* _dst, void* _src, ae_state *_state);
void _mlpcvreport_clear(void* _p);
void _mlpcvreport_destroy(void* _p);
void _smlptrnsession_init(void* _p, ae_state *_state);
void _smlptrnsession_init_copy(void* _dst, void* _src, ae_state *_state);
void _smlptrnsession_clear(void* _p);
void _smlptrnsession_destroy(void* _p);
void _mlpetrnsession_init(void* _p, ae_state *_state);
void _mlpetrnsession_init_copy(void* _dst, void* _src, ae_state *_state);
void _mlpetrnsession_clear(void* _p);
void _mlpetrnsession_destroy(void* _p);
void _mlptrainer_init(void* _p, ae_state *_state);
void _mlptrainer_init_copy(void* _dst, void* _src, ae_state *_state);
void _mlptrainer_clear(void* _p);
void _mlptrainer_destroy(void* _p);
void _mlpparallelizationcv_init(void* _p, ae_state *_state);
void _mlpparallelizationcv_init_copy(void* _dst, void* _src, ae_state *_state);
void _mlpparallelizationcv_clear(void* _p);
void _mlpparallelizationcv_destroy(void* _p);
void pcabuildbasis(/* Real    */ ae_matrix* x,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     /* Real    */ ae_vector* s2,
     /* Real    */ ae_matrix* v,
     ae_state *_state);

}
#endif

