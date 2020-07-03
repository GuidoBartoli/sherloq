# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Copyright (c) 2019 Image Processing Research Group of University Federico II of Naples ('GRIP-UNINA').
# All rights reserved.
# This work should only be used for nonprofit purposes.
#
# By downloading and/or using any of these files, you implicitly agree to all the
# terms of the license, as specified in the document LICENSE.txt
# (included in this package) and online at
# http://www.grip.unina.it/download/LICENSE_OPEN.txt
#
"""
@author: davide.cozzolino
"""

import numpy as np
import numpy.linalg as numpyl
from scipy.ndimage.filters import uniform_filter, maximum_filter
import skimage.morphology as ski

from .feat_spam.spam_np_opt import getSpamRes
from .utility.gaussianMixture import gm

paramSpam_default = {'resTranspose': False, 'uniformQuant': False, \
            'values': np.asarray([-0.8,-0.28,+0.16,+0.7]), \
            'ordCooc': 4, 'mapper': list(), \
            'strides': [8,8], 'numFeat': 4**4, 'radius': [1, 1], \
            'symTranspose': False}
ksize_default = 64
stride_default = 8
satutationProb = 0.95
win_v = 5
win_z = 35

def faetReduce(feat_list, inds, whiteningFlag = False):
    cov_mtx = np.cov(feat_list, rowvar = False, bias = True)
    w, v = np.linalg.eigh(cov_mtx)
    w = w[::-1]
    v = v[:,::-1]
    v = v[:, inds]
    if whiteningFlag:
        v = v / np.sqrt(w[inds])
    return v, w

def getWeights(img, res):
    res_m = uniform_filter(res, (win_v, win_v))
    res_v = uniform_filter(np.square(res), (win_v, win_v)) - np.square(res_m)
    mm = res_v<0.005
    mm[:1,:] = True; mm[-1:,:] = True;
    mm[:,:1] = True; mm[:,-1:] = True;
    
    th_White = 253.0 / 256; rd_White = 3
    sat_mask = ski.binary_opening(img>th_White,ski.disk(rd_White))
    mm = np.logical_or(mm,sat_mask)
    mm = maximum_filter(mm, (win_z,win_z))
    weights = np.logical_not(mm)
    return weights

def getCoocValues(res, img_gray, n_clusters=4, random_state=0):
    #img_gray is the image in range [0,1[
    from sklearn.cluster import KMeans
    weights = getWeights(img_gray, res)
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state).fit(res[weights].reshape((-1,1)))
    values = np.sort(kmeans.cluster_centers_.flatten(), axis=None)
    return values

def getSpamFromNoiseprint(res, img_gray, ksize=ksize_default, stride=stride_default, values=None):
    imgsize = img_gray.shape
    
    paramSpam = dict(paramSpam_default)
    paramSpam['strides'] = (stride,stride)
    if values is not None:
        paramSpam['values'] = values
        
    weights = getWeights(img_gray, res)
    spam, weights, range0, range1 = getSpamRes(res, paramSpam, ksize, weights = weights, paddingModality = 0)
    valid = (weights >= satutationProb)
    
    spam = np.sqrt(np.abs(spam))
    spam   = spam[2:-2,2:-2,:]
    valid  = valid[2:-2,2:-2]
    range0 = range0[2:-2]
    range1 = range1[2:-2]
    
    return spam, valid, range0, range1, imgsize

def EMgu(feats, seed = 0, maxIter = 100, replicates = 10, outliersNlogl = 42):

    randomState = np.random.RandomState(seed)
    gm_data = gm(feats.shape[1], [0,], [2,], outliersProb = 0.01, outliersNlogl = outliersNlogl, dtype=list_valid.dtype)
    gm_data.setRandomParams(feats, regularizer = -1.0, randomState = randomState)
    avrLogl, _, _ = gm_data.EM(feats, maxIter = maxIter, regularizer = -1.0)

    for index in range(1, replicates):

        gm_data_1 = gm(feats.shape[1], [0,], [2,], outliersProb = 0.01, outliersNlogl = outliersNlogl, dtype = list_valid.dtype)
        gm_data_1.setRandomParams(feats, regularizer = -1.0, randomState = randomState)

        avrLogl_1, _, _ = gm_data_1.EM(feats, maxIter = maxIter, regularizer = -1.0)
        if (avrLogl_1>avrLogl):
            gm_data  = gm_data_1
            avrLogl  = avrLogl_1

    _, mahal = gm_data.getNlogl(list_spam)

    mahal = mahal.reshape([shape_spam[0],shape_spam[1],])
    other = dict()
    other['Sigma'] = gm_data.listSigma[0]
    other['mu'] = gm_data.mu
    other['outliersNlogl'] = outliersNlogl
    other['outliersProb'] = gm_data.outliersProb

    return mahal, other

def EMgu_img(spam, valid, extFeat = range(32), seed = 0, maxIter = 100, replicates = 10, outliersNlogl = 42):
    shape_spam = spam.shape
    list_spam  = spam.reshape([shape_spam[0]*shape_spam[1],shape_spam[2]])
    list_valid = list_spam[valid.flatten(),:]
    L, eigs = faetReduce(list_valid, extFeat, True)
    list_spam = np.matmul(list_spam, L)
    list_valid = list_spam[valid.flatten(), :]
    
    randomState = np.random.RandomState(seed)
    gm_data = gm(shape_spam[2], [0,], [2,], outliersProb = 0.01, outliersNlogl = outliersNlogl, dtype=list_valid.dtype)
    gm_data.setRandomParams(list_valid, regularizer = -1.0, randomState = randomState)
    avrLogl, _, _ = gm_data.EM(list_valid, maxIter = maxIter, regularizer = -1.0)

    for index in range(1, replicates):
        gm_data_1 = gm(shape_spam[2], [0,], [2,], outliersProb = 0.01, outliersNlogl = outliersNlogl, dtype = list_valid.dtype)
        gm_data_1.setRandomParams(list_valid, regularizer = -1.0, randomState = randomState)
        avrLogl_1, _, _ = gm_data_1.EM(list_valid, maxIter = maxIter, regularizer = -1.0)
        if (avrLogl_1>avrLogl):
            gm_data  = gm_data_1
            avrLogl  = avrLogl_1
        
    _, mahal = gm_data.getNlogl(list_spam)
    mahal = mahal.reshape([shape_spam[0],shape_spam[1],])
    other = dict()
    other['Sigma'] = gm_data.listSigma[0]
    other['mu'] = gm_data.mu
    other['L'] = L
    other['eigs'] = eigs
    other['outliersNlogl'] = outliersNlogl
    other['outliersProb'] = gm_data.outliersProb
    return mahal, other
