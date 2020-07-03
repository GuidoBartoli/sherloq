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

def getIdemMapper(num):
    return {'num': num,'numIn': num, 'table': np.asarray(range(num),dtype=np.int)}

def getSignSymMapper(occo, n):
    #n = 2 * T + 1
    numIn  = n ** occo
    P = getCombinations(occo, n)
    V = np.ones([numIn, 1], dtype= np.bool)
    table = np.zeros([numIn, 1], dtype= np.int)
    indexOut = 0

    for index in range(numIn):
         if V[index]:
             table[index] = indexOut
             V[index] = False
             H = P[index,:]
             indexS = getPos(n-1-H, n, occo)
             table[indexS] = indexOut
             V[indexS] = False
             H = np.flipud(H)
             indexS = getPos(H, n, occo)
             table[indexS] = indexOut
             V[indexS] = False
             indexS = getPos(n-1-H, n, occo)
             table[indexS] = indexOut
             V[indexS] = False
             indexOut = indexOut + 1

    return {'num': indexOut, 'numIn': numIn, 'table': table}

def getSignMapper(occo, n):
    #n = 2 * T + 1
    numIn  = n ** occo
    numOut = (numIn-1)/2 + 1
    P = getCombinations(occo, n)
    V = np.ones([numIn, 1], dtype= np.bool)
    table = np.zeros([numIn, 1], dtype= np.int)
    indexOut = 0

    for index in range(numIn):
         if V[index]:
             table[index] = indexOut
             V[index] = False
             H = n-1-P[index,:]
             indexS = getPos(H, n, occo)
             table[indexS] = indexOut
             V[indexS] = False
             indexOut = indexOut + 1

    return {'num': numOut, 'numIn': numIn, 'table': table}

def getPos(P, n, occo):
    return np.matmul(P, np.power(n, range(0, occo)))


def getCombinations(occo, n):
    num = n ** occo
    P = np.zeros([num, occo], dtype=np.int)
    P[0, :] = 0
    for indexI in range(1 ,num):
        P[indexI ,:] = P[indexI - 1 ,:]
        for indexJ  in range(occo):
            P[indexI, indexJ] = P[indexI, indexJ] + 1
            if P[indexI, indexJ] >= n:
                P[indexI, indexJ] = 0
            else:
                break
    return P

def mapper2filter(mapper, dtype=np.float32):
    table = mapper['table']
    W = np.zeros([1, 1, len(table), mapper['num']], dtype=dtype)
    for index in range(mapper['num']):
        W[0,0,:,index] = (np.equal(table,index)).squeeze().astype(dtype)
    return W