# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Copyright (c) 2016 Image Processing Research Group of University Federico II of Naples ('GRIP-UNINA').
# All rights reserved.
# This work should only be used for nonprofit purposes.
#
# By downloading and/or using any of these files, you implicitly agree to all the
# terms of the license, as specified in the document LICENSE.txt
# (included in this package) and online at
# http://www.grip.unina.it/download/LICENSE_OPEN.txt
#

import numpy as np
from .residue import getFiltersResidue
from . import mapping as spam_m
from scipy.ndimage.filters import uniform_filter

def quantizerScalarEncoder(x, values):
    y = np.zeros(x.shape, dtype = np.int64)
    th = (values[1:]+values[:-1]) / 2
    for index in range(th.size):
        y += x>th[index]
    return y

def getParams(ordResid, symTranspose, q, T, ordCooc, mapper, strides):
    Wres, Fres, resTranspose = getFiltersResidue(ordResid)
    radius = (np.asarray(Wres.shape[0:2]) - 1) / 2

    n = 2*T + 1
    values = (float(q) * Fres / 256.0) * np.asarray(range(-T,T+1)).astype(np.float)

    radius = radius + (ordCooc - (ordCooc % 2)) / 2
    radius = radius.astype(int)
    
    if isinstance(mapper,dict):
        numFeat = mapper['num']
    elif mapper is 'SignSym':
        mapper = spam_m.getSignSymMapper(ordCooc, n)
        numFeat = mapper['num']
    elif mapper is 'Sign':
        mapper = spam_m.getSignMapper(ordCooc,n)
        numFeat = mapper['num']
    elif mapper is 'Idem':
        mapper = []
        numFeat = (n ** ordCooc)
    else:
        mapper = spam_m.getSignSymMapper(ordCooc, n)
        numFeat = mapper['num']

    if isinstance(strides, (list, tuple)):
        strides = strides[0:2]
    else:
        strides = [strides, strides]

    return {'Wres': Wres, 'resTranspose': resTranspose, 'uniformQuant': True, 'values': values, 'ordCooc': ordCooc, 'mapper': mapper,
            'strides': strides, 'numFeat': numFeat, 'radius': radius,
            'symTranspose': symTranspose}


def computeSpamRes(res, params, weights = list(),  normalize = True):

    ## Quantization & Truncation
    values = params['values']
    resQ = quantizerScalarEncoder(res, values)

    ## Coocorance
    ordCooc = params['ordCooc']
    n = (params['values']).size
    dim = int(ordCooc + 1 - (ordCooc % 2))
    indexL = int((dim - 1) / 2)

    shapeR = np.asarray(resQ.shape[:2]) - dim + 1
    resH = np.zeros(shapeR, dtype = np.int)
    resV = np.zeros(shapeR, dtype = np.int)

    for indexP in range(ordCooc):
        nn = (n ** indexP)
        resH += resQ[indexL:(shapeR[0] + indexL), indexP:(shapeR[1] + indexP)] * nn
        resV += resQ[indexP:(shapeR[0] + indexP), indexL:(shapeR[1] + indexL)] * nn

    ## Mappeing
    mapper = params['mapper']
    if len(mapper) > 0:
        resH = mapper['table'][resH].squeeze()
        resV = mapper['table'][resV].squeeze()

    ## Hist
    strides = params['strides']
    numFeat = max(max(np.max(resH), np.max(resV))+1,params['numFeat'])

    shapeR = resH.shape
    range0 = np.arange(0, shapeR[0]-strides[0]+1, strides[0], dtype=np.uint16)
    range1 = np.arange(0, shapeR[1]-strides[1]+1, strides[1], dtype=np.uint16)
    rangeH = np.arange(0, numFeat+2, dtype = resH.dtype) # range(0, numFeat+1)
    if normalize:
        out_dtype = np.float32
    else:
        out_dtype = np.uint32
        
    spamH = np.zeros([range0.size, range1.size, numFeat+1], dtype = out_dtype)
    spamV = np.zeros([range0.size, range1.size, numFeat+1], dtype = out_dtype)

    if len(weights) > 0:
        weights = weights[indexL:(shapeR[0] + indexL), indexL:(shapeR[1] + indexL)] ## clip weights
        resH[np.logical_not(weights)] = numFeat
        resV[np.logical_not(weights)] = numFeat
        weights = weights.astype(dtype=out_dtype)
    else:
        weights = np.ones(resH.shape, dtype=out_dtype)
        
    for index0 in range(range0.size):
        for index1 in range(range1.size):
            pos0 = range0[index0]
            end0 = strides[0] + pos0
            pos1 = range1[index1]
            end1 = strides[1] + pos1

            spamH[index0, index1, :], _ = np.histogram(resH[pos0:end0, pos1:end1], rangeH, density=False)
            spamV[index0, index1, :], _ = np.histogram(resV[pos0:end0, pos1:end1], rangeH, density=False)

    spamW = (strides[0]*strides[1]) - spamH[:,:,-1]
    spamH = spamH[:,:,:-1]
    spamV = spamV[:,:,:-1]
    if normalize:
        spamH = spamH / np.maximum(spamW[:,:,np.newaxis],1e-20)
        spamV = spamV / np.maximum(spamW[:,:,np.newaxis],1e-20)
        spamW = spamW / (strides[0]*strides[1])
        
    spam = np.concatenate([spamH, spamV], 2)

    ## Border
    range0 = range0 + indexL + (strides[0] - 1.0) / 2.0
    range1 = range1 + indexL + (strides[1] - 1.0) / 2.0

    return spam, spamW, range0, range1


def getSpamRes(res, params, ksize, weights = list(), paddingModality = 0):

    strides = params['strides']
    if isinstance(ksize, (list,tuple)):
        ksize = ksize[0:2]
    else:
        ksize = [ksize, ksize]

    ksize[0] = int(ksize[0] / strides[0])
    ksize[1] = int(ksize[1] / strides[1])
    spam, spamW, range0, range1 = computeSpamRes(res, params, weights = weights, normalize=True)
    
    spamW = np.maximum(uniform_filter(spamW, (ksize[0], ksize[1],  ), mode='constant', cval=0.0), 0.0)
    spam  = np.maximum(uniform_filter(spam , (ksize[0], ksize[1], 1), mode='constant', cval=0.0), 0.0)
    spam  = spam / np.maximum(spamW[:,:,np.newaxis], 1e-20)
    
    if paddingModality == 0:
        ind0f = int(np.floor((ksize[0] - 1.0) / 2.0))
        ind0c = int(np.ceil( (ksize[0] - 1.0) / 2.0))
        ind1f = int(np.floor((ksize[1] - 1.0) / 2.0))
        ind1c = int(np.ceil( (ksize[1] - 1.0) / 2.0))
        range0 = (range0[ind0f:-ind0c] + range0[ind0c:-ind0f])/2.0
        range1 = (range1[ind1f:-ind1c] + range1[ind1c:-ind1f])/2.0
        spamW  = spamW[ind0c:-ind0f,ind1c:-ind1f]
        spam   = spam [ind0c:-ind0f,ind1c:-ind1f,:]
    
    return spam, spamW, range0, range1
