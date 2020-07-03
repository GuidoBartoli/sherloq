# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# Copyright (c) 2018 Image Processing Research Group of University Federico II of Naples ('GRIP-UNINA').
# All rights reserved.
# This work should only be used for nonprofit purposes.
#
# By downloading and/or using any of these files, you implicitly agree to all the
# terms of the license, as specified in the document LICENSE.txt
# (included in this package) and online at
# http://www.grip.unina.it/download/LICENSE_OPEN.txt
#

from PIL.JpegImagePlugin import convert_dict_qtables
from PIL import Image
import numpy as np


def imread2f_pil(stream, channel = 1, dtype = np.float32):
    img = Image.open(stream)
    mode = img.mode
    
    if channel == 3:
        img = img.convert('RGB')
        img = np.asarray(img).astype(dtype) / 256.0
    elif channel == 1:
        if img.mode == 'L':
            img = np.asarray(img).astype(dtype) / 256.0
        else:
            img = img.convert('RGB')
            img = np.asarray(img).astype(dtype)
            img = (0.299 * img[:, :, 0] + 0.587 * img[:, :, 1] + 0.114 * img[:, :, 2])/256.0
    else:
        img = np.asarray(img).astype(dtype) / 256.0
    return img, mode
try:
    import rawpy
    def imread2f_raw(stream, channel = 1, dtype = np.float32):
        raw = rawpy.imread(stream)
        img = raw.postprocess()
        raw.close()
        ori_dtype = img.dtype
        img = np.asarray(img).astype(dtype)
        if channel == 1:
            img = (0.299 * img[:, :, 0] + 0.587 * img[:, :, 1] + 0.114 * img[:, :, 2])

        if ori_dtype == np.uint8:
            img = img / 256.0
        elif ori_dtype == np.uint16:
            img = img / (2.0 ** 16)
        elif ori_dtype == np.uint32:
            img = img / (2.0 ** 32)
        elif ori_dtype == np.uint64:
            img = img / (2.0 ** 64)
        elif ori_dtype == np.uint128:
            img = img / (2.0 ** 128)
        return img, 'RAW'
except:
    pass

def imread2f(stream, channel = 1, dtype = np.float32):
    try:
       return imread2f_raw(stream, channel=channel, dtype=dtype)
    except:
       return imread2f_pil(stream, channel=channel, dtype=dtype)
    
    
def jpeg_qtableinv(stream, tnum=0, force_baseline=None):
    assert tnum == 0 or tnum == 1, 'Table number must be 0 or 1'

    if force_baseline is None:
        th_high = 32767
    elif force_baseline == 0:
        th_high = 32767
    else:
        th_high = 255

    h = np.asarray(convert_dict_qtables(Image.open(stream).quantization)[tnum]).reshape((8, 8))

    if tnum == 0:
        # This is table 0 (the luminance table):
        t = np.matrix(
             [[16,  11,  10,  16,  24,  40,  51,  61],
              [12,  12,  14,  19,  26,  58,  60,  55],
              [14,  13,  16,  24,  40,  57,  69,  56],
              [14,  17,  22,  29,  51,  87,  80,  62],
              [18,  22,  37,  56,  68, 109, 103,  77],
              [24,  35,  55,  64,  81, 104, 113,  92],
              [49,  64,  78,  87, 103, 121, 120, 101],
              [72,  92,  95,  98, 112, 100, 103,  99]])

    elif tnum == 1:
        # This is table 1 (the chrominance table):
        t = np.matrix(
            [[17,  18,  24,  47,  99,  99,  99,  99],
             [18,  21,  26,  66,  99,  99,  99,  99],
             [24,  26,  56,  99,  99,  99,  99,  99],
             [47,  66,  99,  99,  99,  99,  99,  99],
             [99,  99,  99,  99,  99,  99,  99,  99],
             [99,  99,  99,  99,  99,  99,  99,  99],
             [99,  99,  99,  99,  99,  99,  99,  99],
             [99,  99,  99,  99,  99,  99,  99,  99]])

    else:
        raise ValueError(tnum, 'Table number must be 0 or 1')

    h_down = np.divide((2 * h-1), (2 * t))
    h_up   = np.divide((2 * h+1), (2 * t))
    if np.all(h == 1): return 100
    x_down = (h_down[h > 1]).max()
    x_up   = (h_up[h < th_high]).min() if (h < th_high).any() else None
    if x_up is None:
        s = 1
    elif x_down > 1 and x_up > 1:
        s = np.ceil(50 / x_up)
    elif x_up < 1:
        s = np.ceil(50*(2 - x_up))
    else:
        s = 50
    return s


from scipy.interpolate import interp1d
def resizeMapWithPadding(x, range0, range1, shapeOut):
    range0 = range0.flatten()
    range1 = range1.flatten()
    xv = np.arange(shapeOut[1])
    yv = np.arange(shapeOut[0])
    y = interp1d(range1, x    , axis=1, kind='nearest', fill_value='extrapolate', assume_sorted=True, bounds_error=False)
    y = interp1d(range0, y(xv), axis=0, kind='nearest', fill_value='extrapolate', assume_sorted=True, bounds_error=False)
    return y(yv).astype(x.dtype)


def computeMetricsContinue(values, gt0, gt1):
    values = values.flatten().astype(np.float32)
    gt0 = gt0.flatten().astype(np.float32)
    gt1 = gt1.flatten().astype(np.float32)
    inds = np.argsort(values)
    inds = inds[(gt0[inds]+gt1[inds])>0]
    vet_th = values[inds]
    gt0 = gt0[inds]
    gt1 = gt1[inds]
    
    TN = np.cumsum(gt0)
    FN = np.cumsum(gt1)
    FP = np.sum(gt0) - TN
    TP = np.sum(gt1) - FN

    return FP, TP, FN, TN, vet_th 

def computeMCC(values, gt0, gt1):
    FP, TP, FN, TN, vet_th  = computeMetricsContinue(values, gt0, gt1)
    mcc = np.abs(TP*TN - FP*FN) / np.maximum(np.sqrt((TP + FP)*(TP + FN)*(TN+ FP)*(TN + FN) ), 1e-32)
    return mcc, vet_th