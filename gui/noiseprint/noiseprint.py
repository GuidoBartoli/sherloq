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
import tensorflow.compat.v1 as tf

tf.disable_v2_behavior()
import os.path
from .network import FullConvNet

slide = 1024  # 3072
largeLimit = 1050000  # 9437184
overlap = 34

chkpt_folder = os.path.join(os.path.dirname(__file__), "./nets/%s_jpg%d/model")
tf.reset_default_graph()
x_data = tf.placeholder(tf.float32, [1, None, None, 1], name="x_data")
net = FullConvNet(x_data, 0.9, tf.constant(False), num_levels=17)
saver = tf.train.Saver(net.variables_list)

configSess = tf.ConfigProto()
configSess.gpu_options.allow_growth = True
# configSess = tf.ConfigProto(gpu_options=tf.GPUOptions(per_process_gpu_memory_fraction=0.95))


def genNoiseprint(img, QF=101, model_name="net"):
    if QF > 100:
        QF = 101
    chkpt_fname = chkpt_folder % (model_name, QF)

    with tf.Session(config=configSess) as sess:
        saver.restore(sess, chkpt_fname)

        if img.shape[0] * img.shape[1] > largeLimit:
            # print(' %dx%d large %3d' % (img.shape[0], img.shape[1], QF))
            # for large image the network is executed windows with partial overlapping
            res = np.zeros((img.shape[0], img.shape[1]), np.float32)
            for index0 in range(0, img.shape[0], slide):
                index0start = index0 - overlap
                index0end = index0 + slide + overlap

                for index1 in range(0, img.shape[1], slide):
                    index1start = index1 - overlap
                    index1end = index1 + slide + overlap
                    clip = img[
                        max(index0start, 0) : min(index0end, img.shape[0]),
                        max(index1start, 0) : min(index1end, img.shape[1]),
                    ]
                    resB = sess.run(net.output, feed_dict={x_data: clip[np.newaxis, :, :, np.newaxis]})
                    resB = np.squeeze(resB)

                    if index0 > 0:
                        resB = resB[overlap:, :]
                    if index1 > 0:
                        resB = resB[:, overlap:]
                    resB = resB[: min(slide, resB.shape[0]), : min(slide, resB.shape[1])]

                    res[index0 : min(index0 + slide, res.shape[0]), index1 : min(index1 + slide, res.shape[1])] = resB
        else:
            # print(' %dx%d small %3d' % (img.shape[0], img.shape[1], QF))
            res = sess.run(net.output, feed_dict={x_data: img[np.newaxis, :, :, np.newaxis]})
            res = np.squeeze(res)
    return res
