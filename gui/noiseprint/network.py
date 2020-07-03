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

import numpy as np
import tensorflow as tf

class FullConvNet(object):
  """FullConvNet model."""
  
  def __init__(self, images, bnorm_decay, falg_train, num_levels = 17, padding = 'SAME'):
    """FullConvNet constructor."""
    
    self._num_levels = num_levels
    self._actfun   = [tf.nn.relu, ] * (self._num_levels-1) + [tf.identity, ]
    self._f_size   = [3, ] * self._num_levels
    self._f_num    = [64, ] *(self._num_levels-1) + [1, ]
    self._f_stride = [1, ] * self._num_levels
    self._bnorm    = [False, ] + [True,]*(self._num_levels-2) + [False,]
    self._res      = [0, ] * self._num_levels
    self._bnorm_init_var = 1e-4
    self._bnorm_init_gamma = np.sqrt(2.0/(9.0*64.0))
    self._bnorm_epsilon = 1e-5
    self._bnorm_decay = bnorm_decay
    
    self.level = [None, ] * self._num_levels
    self.input = images
    self.falg_train = falg_train
    self.extra_train = []
    self.variables_list = []
    self.trainable_list = []
    self.decay_list = []
    self.padding = padding
    
    x = self.input
    for i in range(self._num_levels):
      with tf.variable_scope('level_%d' % i):
        x = self._conv(x, self._f_size[i], self._f_num[i], self._f_stride[i], name = 'conv')      
        if self._bnorm[i]:
            x = self._batch_norm(x, name = 'bn')        
        x = self._bias(x, name = 'bias')
        if self._res[i]>0:
            x = x + self.level[i-self._res[i]]          
        x = self._actfun[i](x, name = 'active')
        self.level[i] = x     
    self.output = x
  
  def _batch_norm(self, x, name = 'bnorm'):
    """Batch normalization."""
    with tf.variable_scope(name):
      params_shape = [x.get_shape()[-1]]
      
      moving_mean = tf.get_variable(
            'moving_mean', params_shape, tf.float32,
            initializer=tf.constant_initializer(0.0, dtype=tf.float32),
            trainable=False)
      moving_variance = tf.get_variable(
            'moving_variance', params_shape, tf.float32,
            initializer=tf.constant_initializer(self._bnorm_init_var, dtype=tf.float32),
            trainable=False)
      self.variables_list.append(moving_mean)
      self.variables_list.append(moving_variance)
      
      gamma = tf.get_variable(
          'gamma', params_shape, tf.float32,
          initializer=tf.random_normal_initializer(stddev=self._bnorm_init_gamma, dtype=tf.float32))
      self.variables_list.append(gamma)
      self.trainable_list.append(gamma)
 
      local_mean, local_variance = tf.nn.moments(x, [0, 1, 2], name='moments')
      
      mean, variance = tf.cond(
        self.falg_train, lambda: (local_mean, local_variance),
        lambda: (moving_mean, moving_variance))
        
      self.extra_train.append(moving_mean.assign_sub((1.0 - self._bnorm_decay) * (moving_mean - local_mean)))
      self.extra_train.append(moving_variance.assign_sub((1.0 - self._bnorm_decay) * (moving_variance - local_variance)))
         
      y = tf.nn.batch_normalization(
          x, mean, variance, None, gamma, self._bnorm_epsilon)
      y.set_shape(x.get_shape())
    return y
  
  def _bias(self, x, name = 'bias'):
    """Bias term."""
    with tf.variable_scope(name):
      params_shape = [x.get_shape()[-1]]
      beta = tf.get_variable(
          'beta', params_shape, tf.float32,
          initializer=tf.constant_initializer(0.0, dtype=tf.float32))
      self.variables_list.append(beta)
      self.trainable_list.append(beta)
      y = x + beta
    return y 
  
  def _conv(self, x, filter_size, out_filters, stride, name='conv'):
    """Convolution."""
    with tf.variable_scope(name):
      in_filters = int(x.get_shape()[-1])
      n = filter_size * filter_size * np.maximum(in_filters, out_filters)
      kernel = tf.get_variable(
          'weights', [filter_size, filter_size, in_filters, out_filters],
          tf.float32, initializer=tf.random_normal_initializer(
              stddev=np.sqrt(2.0/n), dtype=tf.float32))
      self.variables_list.append(kernel)
      self.trainable_list.append(kernel)
      self.decay_list.append(kernel)
      y = tf.nn.conv2d(x, kernel, [1, stride, stride, 1], padding=self.padding)
    return y
