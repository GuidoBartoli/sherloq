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
from .noiseprint import genNoiseprint
from .post_em import EMgu_img, getSpamFromNoiseprint
from .utility.utilityRead import resizeMapWithPadding
from .utility.utilityRead import imread2f
from .utility.utilityRead import jpeg_qtableinv


def noiseprint_blind_file(filename, model_name="net"):
    try:
        img, mode = imread2f(filename, channel=1)
    except:
        # print('Error opening image')
        return -1, -1, -1e10, None, None, None, None, None, None

    try:
        QF = jpeg_qtableinv(filename)
        # print('QF=', QF)
    except:
        QF = 200

    mapp, valid, range0, range1, imgsize, other = noiseprint_blind(img, QF, model_name=model_name)
    return QF, mapp, valid, range0, range1, imgsize, other


def noiseprint_blind(img, QF, model_name="net"):
    res = genNoiseprint(img, QF, model_name)
    assert img.shape == res.shape
    return noiseprint_blind_post(res, img)


def noiseprint_blind_post(res, img):
    spam, valid, range0, range1, imgsize = getSpamFromNoiseprint(res, img)

    if np.sum(valid) < 50:
        # print('error too small %d' % np.sum(weights))
        return None, valid, range0, range1, imgsize, dict()

    mapp, other = EMgu_img(spam, valid, extFeat=range(32), seed=0, maxIter=100, replicates=10, outliersNlogl=42)

    return mapp, valid, range0, range1, imgsize, other


def genMappFloat(mapp, valid, range0, range1, imgsize):
    mapp_s = np.copy(mapp)
    mapp_s[valid == 0] = np.min(mapp_s[valid > 0])

    mapp_s = resizeMapWithPadding(mapp_s, range0, range1, imgsize)

    return mapp_s


def genMappUint8(mapp, valid, range0, range1, imgsize, vmax=None, vmin=None):
    mapp_s = np.copy(mapp)
    mapp_s[valid == 0] = np.min(mapp_s[valid > 0])

    if vmax is None:
        vmax = np.nanmax(mapp_s)
    if vmin is None:
        vmin = np.nanmin(mapp_s)

    mapUint8 = (255 * (mapp_s.clip(vmin, vmax) - vmin) / (vmax - vmin)).clip(0, 255).astype(np.uint8)
    mapUint8 = 255 - resizeMapWithPadding(mapUint8, range0, range1, imgsize)

    return mapUint8
