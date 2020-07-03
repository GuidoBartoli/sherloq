# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Copyright (c) 2017 Image Processing Research Group of University Federico II of Naples ('GRIP-UNINA').
# This software is delivered with Government Purpose Rights (GPR) under agreement number FA8750-16-2-0204.
#
# By downloading and/or using any of these files, you implicitly agree to all the
# terms of the license, as specified in the document LICENSE.txt
# (included in this package) 
#

import numpy as np
from scipy.linalg import eigvalsh
from numpy.linalg import cholesky
from numpy.linalg import eigh

class gm:
    prioriProb = 0
    outliersProb = 0
    outliersNlogl = 0
    mu = 0
    listSigma = []
    listSigmaInds = []
    listSigmaType = []

    # sigmaType = 0 # isotropic covariance
    # sigmaType = 1 # diagonal covariance
    # sigmaType = 2 # full covariance

    # outliersProb <  0 # outliers are not managed
    # outliersProb >= 0 # outliers are managed throught fixed nlogl (negative log likelihood)
    # TODO: outliers managed throught fixed probability

    def __init__(self, dim, listSigmaInds, listSigmaType, outliersProb = -1, outliersNlogl = 0, dtype = np.float32):
        K = len(listSigmaInds)
        S = len(listSigmaType)

        self.listSigmaInds = listSigmaInds
        self.listSigmaType = listSigmaType
        self.outliersProb  = outliersProb
        self.outliersNlogl = outliersNlogl
        self.prioriProb = (1.0-self.outliersProb) * np.ones((K, 1), dtype=dtype) / K
        self.mu = np.zeros((K, dim), dtype=dtype)
        self.listSigma = [None, ] * S

        for s in range(S):
            sigmaType = self.listSigmaType[s]
            if sigmaType == 2:  # full covariance
                self.listSigma[s] = np.ones([dim, dim], dtype = dtype)
            elif sigmaType == 1:  # diagonal covariance
                self.listSigma[s] = np.ones([1, dim], dtype = dtype)
            else:
                self.listSigma[s] = np.ones([], dtype = dtype)

    def setRandomParams(self, X, regularizer = 0, randomState = np.random.get_state()):
        [N, dim] = X.shape
        K = len(self.listSigmaInds)
        S = len(self.listSigmaType)
        dtype = X.dtype

        if self.outliersProb > 0:
            self.prioriProb = (1.0-self.outliersProb) * np.ones((K, 1), dtype=dtype) / K
        else:
            self.prioriProb = np.ones((K, 1), dtype=dtype) / K

        inds = randomState.random_integers(low=0,high=(N-1),size=(K,))
        self.mu = X[inds, :]
        varX = np.var(X, axis=0, keepdims=True)
        if regularizer>0:
            varX = varX + regularizer
        elif regularizer<0:
            varX = varX  + np.abs(regularizer*np.spacing(np.max(varX)))

        for s in range(S):
            sigmaType = self.listSigmaType[s]
            if sigmaType == 2:  # full covariance
                self.listSigma[s] = np.diag(varX.flatten())
            elif sigmaType == 1:  # diagonal covariance
                self.listSigma[s] = varX
            else:
                self.listSigma[s] = np.mean(varX)
        return inds

    def setRandomParamsW(self, X, weights, regularizer = 0, randomState = np.random.get_state(), meanFlag = False):
        [N, dim] = X.shape
        K = len(self.listSigmaInds)
        S = len(self.listSigmaType)
        dtype = X.dtype

        if self.outliersProb > 0:
            self.prioriProb = (1.0-self.outliersProb) * np.ones((K, 1), dtype=dtype) / K
        else:
            self.prioriProb = np.ones((K, 1), dtype=dtype) / K

        avrX = np.mean(X*weights, axis=0, keepdims=True)/np.mean(weights)
        varX = np.mean(weights *((X - avrX) ** 2), axis=0, keepdims=True)/np.mean(weights)

        indsW = np.sum(weights)*randomState.random_sample(size=(K,))
        inds = [None, ] * K
        weights = np.cumsum(weights.flatten())
        for index in range(K):
            inds[index] = np.count_nonzero(weights<=indsW[index])

        self.mu = X[inds, :]
        if meanFlag: self.mu[0,:] = avrX
        #varX = np.var(X, axis=0, keepdims=True)
        if regularizer>0:
            varX = varX + regularizer
        elif regularizer<0:
            varX = varX + np.abs(regularizer*np.spacing(np.max(varX)))

        for s in range(S):
            sigmaType = self.listSigmaType[s]
            if sigmaType == 2:  # full covariance
                self.listSigma[s] = np.diag(varX.flatten())
            elif sigmaType == 1:  # diagonal covariance
                self.listSigma[s] = varX
            else:
                self.listSigma[s] = np.mean(varX)
        return inds

    def getNlogl(self, X):
        [N, dim] = X.shape
        K = len(self.listSigmaInds)
        S = len(self.listSigmaType)
        dtype = X.dtype

        K0 = K
        if self.outliersProb >= 0: K0 = K+1

        nlogl = np.zeros([N, K0], dtype = dtype)
        mahal = np.zeros([N, K ], dtype = dtype)
        listLogDet = [None, ] * S
        listLowMtx = [None, ] * S
        for s in range(S):
            sigmaType = self.listSigmaType[s]
            sigma = self.listSigma[s]
            if sigmaType == 2:  # full covariance
                try:
                    listLowMtx[s] = cholesky(sigma)
                except:
                    # exceptional regularization
                    sigma_w, sigma_v = eigh(np.real(sigma))
                    sigma_w = np.maximum(sigma_w, np.spacing(np.max(sigma_w)))
                    sigma = np.matmul(np.matmul(sigma_v, np.diag(sigma_w)), (np.transpose(sigma_v,[1,0])))
                    try:
                        listLowMtx[s] = cholesky(sigma)
                    except:
                        sigma_w, sigma_v = eigh(np.real(sigma))
                        sigma_w = np.maximum(sigma_w, np.spacing(np.max(sigma_w)))
                        #print(np.min(sigma_w))
                        sigma = np.matmul(np.matmul(sigma_v, np.diag(sigma_w)), (np.transpose(sigma_v,[1,0])))
                        #print(sigma)
                        listLowMtx[s] = cholesky(sigma)
                diagLowMtx = np.diag(listLowMtx[s])
                listLogDet[s] = 2 * np.sum(np.log(diagLowMtx))
            elif sigmaType == 1:  # diagonal covariance
                listLowMtx[s] = np.sqrt(sigma)
                listLogDet[s] = np.sum(np.log(sigma))
            else: # isotropic covariance
                listLowMtx[s] = np.sqrt(sigma)
                listLogDet[s] = dim * np.log(sigma)

        constPi = dim*np.log(2*np.pi)
        for k in range(K):
            s = self.listSigmaInds[k]
            sigmaType = self.listSigmaType[s]
            lowMtx = listLowMtx[s]
            logDet = listLogDet[s]

            Xmu =  X - self.mu[k,:]

            if sigmaType == 2:  # full covariance
                Xmu = np.linalg.solve(lowMtx, Xmu.transpose()).transpose()
            elif sigmaType == 1:  # diagonal covariance
                Xmu = Xmu / lowMtx
            else:  # isotropic covariance
                Xmu = Xmu / lowMtx

            mahal[:,k] = np.sum(Xmu * Xmu, axis = 1)

            nlogl[:,k] = 0.5 * (mahal[:,k] + logDet + constPi)

        if self.outliersProb >= 0:
            nlogl[:, K] = self.outliersNlogl

        return nlogl, mahal

    def getLoglh(self, X):
        nlogl, _ = self.getNlogl(X)
        logPrb = np.log(self.prioriProb)
        if self.outliersProb >= 0:
            #print(self.outliersProb)
            logPrb = np.append(logPrb.squeeze(), np.log(self.outliersProb))
            logPrb = logPrb.reshape((-1,1))

        return logPrb.transpose((1,0)) - nlogl

    def getLoglhInlier(self, X):
        nlogl, _ = self.getNlogl(X)
        K = self.prioriProb.size
        logPrb = np.log(self.prioriProb)
        logit = logPrb.transpose((1, 0)) - nlogl[:, :K]

        maxll = np.max(logit, axis=1, keepdims=True)
        prob = np.exp(logit - maxll)
        dem = np.sum(prob, axis=1, keepdims=True)
        #return (np.log(dem) + maxll - np.log(np.sum(self.prioriProb)))
        return (np.log(dem) + maxll - np.log(np.sum(self.outliersProb)))

    def maximizationParam(self, X, post, regularizer = 0):
        [N, dim] = X.shape
        K = len(self.listSigmaInds)
        S = len(self.listSigmaType)
        dtype = X.dtype

        self.prioriProb = np.sum(post[:,:K], axis=0, keepdims=True).transpose([1, 0])

        self.mu = np.tensordot(post, X, (0, 0)) / self.prioriProb
        for s in range(S):
            sigmaType = self.listSigmaType[s]
            if sigmaType == 2:  # full covariance
                sigma = np.zeros([dim, dim], dtype=dtype)
                sigmadem = np.zeros([], dtype=dtype)
                for k in range(K):
                    if s == self.listSigmaInds[k]:
                        Xmu = X - self.mu[(k,), :]
                        Xmu = np.sqrt(post[:, (k,)]) * Xmu
                        sigma = sigma + np.tensordot(Xmu, Xmu, (0, 0))
                        sigmadem += self.prioriProb[k, 0]
                sigma = sigma / sigmadem
                if regularizer > 0:
                    sigma = sigma + regularizer * np.eye(dim)
                elif regularizer < 0:
                    #sigma = sigma - regularizer * np.spacing(np.max(np.linalg.eigvalsh(sigma))) * np.eye(dim)
                    sigma = sigma + np.abs(regularizer * np.spacing(eigvalsh(sigma, eigvals=(dim - 1, dim - 1)))) * np.eye(dim)
            elif sigmaType == 1:  # diagonal covariance
                sigma = np.zeros([1, dim], dtype=dtype)
                sigmadem = np.zeros([], dtype=dtype)
                for k in range(K):
                    if s == self.listSigmaInds[k]:
                        Xmu = X - self.mu[(k,), :]
                        sigma = sigma + np.tensordot(post[:, (k,)], (Xmu * Xmu), (0, 0))
                        sigmadem += self.prioriProb[k, 0]
                sigma = sigma / sigmadem
                if regularizer > 0:
                    sigma = sigma + regularizer
                elif regularizer < 0:
                    sigma = sigma + + np.abs(regularizer * np.spacing(np.max(sigma)))
            else:  # isotropic covariance
                sigma = np.zeros([], dtype=dtype)
                sigmadem = np.zeros([], dtype=dtype)
                for k in range(K):
                    if s == self.listSigmaInds[k]:
                        Xmu = X - self.mu[(k,), :]
                        sigma = sigma + np.dot(post[:, k], np.mean((Xmu * Xmu), axis=1))
                        sigmadem += self.prioriProb[k, 0]
                sigma = sigma / sigmadem
                if regularizer > 0:
                    sigma = sigma + regularizer
                elif regularizer < 0:
                    sigma = sigma + np.abs(regularizer * np.spacing(sigma))
            self.listSigma[s] = sigma

        # normalize PComponents
        if self.outliersProb < 0:
            self.prioriProb = self.prioriProb / np.sum(self.prioriProb)
        else:
            self.outliersProb = np.sum(post[:,K])
            dem = self.outliersProb +  np.sum(self.prioriProb)
            self.prioriProb = self.prioriProb / dem
            self.outliersProb = self.outliersProb / dem

    def expectation(self, X):
        [post, avrLogl] = softmax(self.getLoglh(X))
        return post, avrLogl

    def expectationWeighed(self, X, weighed):
        [post, avrLogl] = softmaxWeighed(self.getLoglh(X), weighed)
        return post, avrLogl

    def MEstep(self, X, post, regularizer = 0):
        self.maximizationParam(X, post, regularizer = regularizer)
        [post, avrLogl] = self.expectation(X)
        return post, avrLogl

    def MEstepWeighed(self, X, weights, post, regularizer = 0):
        self.maximizationParam(X, post * weights, regularizer = regularizer)
        [post, avrLogl] = self.expectationWeighed(X, weights)
        return post, avrLogl

    def EM(self, X, regularizer, maxIter, relErr = 1e-5):
        [post, avrLogl_old] = self.expectation(X)

        flagExit = 1
        # flagExit = 1 # max number of iteretions
        # flagExit = 0 # converged
        for iter in range(maxIter):
            [post, avrLogl] = self.MEstep(X, post, regularizer = regularizer)

            diff = avrLogl - avrLogl_old
            if (diff >= 0) & (diff < relErr * np.abs(avrLogl)):
                flagExit = 0
                break
            avrLogl_old = avrLogl

        return avrLogl, flagExit, iter


    def EMweighed(self, X, weights, regularizer, maxIter, relErr=1e-5):
        [post, avrLogl_old] = self.expectationWeighed(X, weights)

        flagExit = 1
        # flagExit = 1 # max number of iteretions
        # flagExit = 0 # converged
        for iter in range(maxIter):
            [post, avrLogl] = self.MEstepWeighed(X, weights, post, regularizer=regularizer)

            diff = avrLogl - avrLogl_old
            if (diff >= 0) & (diff < relErr * np.abs(avrLogl)):
                flagExit = 0
                break
            avrLogl_old = avrLogl

        return avrLogl, flagExit, iter

def softmax(logit):
    maxll = np.max(logit, axis = 1, keepdims=True)
    prob  = np.exp(logit - maxll)
    dem   = np.sum(prob, axis = 1, keepdims=True)
    prob  = prob / dem
    avrLogl = np.mean(np.log(dem) + maxll)
    return prob, avrLogl

def softmaxWeighed(logit, weights):
    maxll = np.max(logit, axis = 1, keepdims=True)
    prob  = np.exp(logit - maxll)
    dem   = np.sum(prob, axis = 1, keepdims=True)
    prob  = prob / dem
    avrLogl = np.mean(weights * (np.log(dem) + maxll)) /  np.mean(weights)
    return prob, avrLogl
