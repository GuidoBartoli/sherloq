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
from .mapping import getCombinations

def getFiltersResidue(res):
    try:
       if not(isinstance(res,basestring)): res = str(res)
    except: 
       if not(isinstance(res,str)): res = str(res)
    resTranspose = False

    if res == '0':
        # zero order
        W = np.zeros([3,3,1,1],dtype=np.float32)
        W[1, 1, 0, 0] = 1.0
        F = 1.0
    elif res == '1':
        # 1st order
        W = np.zeros([3,3,1,1],dtype=np.float32)
        W[1, :, 0, 0] = [-1.0, 1.0, 0]
        resTranspose = True
        F = 1.0
    elif res == '2':
        # 2nd order
        W = np.zeros([3, 3, 1, 1], dtype=np.float32)
        W[1, :, 0, 0] = [-1.0, 2.0, -1.0]
        resTranspose = True
        F = 2.0
    elif res == '3':
        # 3rd order
        W = np.zeros([5, 5, 1, 1], dtype=np.float32)
        W[2, :, 0, 0] = [ 0.0, +1.0, -3.0, +3.0, -1.0]
        resTranspose = True
        F = 3.0
    elif res == '5x5':
        # 5x5
        W = np.zeros([5, 5, 1, 1], dtype=np.float32)
        W[0, :, 0, 0] = [-1.0, +2.0,  -2.0, +2.0, -1.0]
        W[1, :, 0, 0] = [+2.0, -6.0,  +8.0, -6.0, +2.0]
        W[2, :, 0, 0] = [-2.0, +8.0, -12.0, +8.0, -2.0]
        W[3, :, 0, 0] = [+2.0, -6.0,  +8.0, -6.0, +2.0]
        W[4, :, 0, 0] = [-1.0, +2.0,  -2.0, +2.0, -1.0]
        F = 12.0
    else:
        W = np.zeros([1, 1, 1, 1], dtype=np.float32)
        F = 0

    # W = [1, -4, 6, -4, 1;] / 6;
    # W = [- 1, 5, -10, 10, -5, 1;] / 10;
    # W = [- 1, 6, -15, 20, -15, 6, -1;] / 20;

    return W, F, resTranspose

def getFilterOcco(occo, values):
    dim = int(occo + 1 - (occo%2))
    values = np.asarray(values).astype(np.float32)
    M = np.zeros([occo, dim * dim], dtype=np.float32)
    f = int((dim-1)/2)

    for index in range(occo):
        R = np.zeros([dim, dim], dtype=np.float32)
        R[index, f] = 1.0
        M[index, :] = R.flatten()

    P = getCombinations(occo, len(values))
    V = values[P]
    n = V.shape[0]

    W =  np.matmul(V, M)
    B = -0.5 * np.sum(np.square(V), axis=1)

    W = np.reshape(W, [n, dim, dim, 1])
    W = np.transpose(W, [1, 2, 3, 0])

    return W, B
