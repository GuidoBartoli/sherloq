from time import time

import cv2 as cv
import numpy as np
from PySide2.QtWidgets import (
    QProgressDialog,
    QVBoxLayout,
    QHBoxLayout,
    QComboBox,
    QSpinBox,
    QLabel)

from tools import ToolWidget
from utility import elapsed_time, normalize_mat
from viewer import ImageViewer


class MinMaxWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(MinMaxWidget, self).__init__(parent)

        self.chan_combo = QComboBox()
        self.chan_combo.addItems(
            [self.tr('Luminance'), self.tr('Red'), self.tr('Green'), self.tr('Blue'), self.tr('RGB Norm')])

        colors = [self.tr('Red'), self.tr('Green'), self.tr('Blue'), self.tr('White'), self.tr('Black')]

        self.min_combo = QComboBox()
        self.min_combo.addItems(colors)
        self.min_combo.setCurrentIndex(1)

        self.max_combo = QComboBox()
        self.max_combo.addItems(colors)
        self.max_combo.setCurrentIndex(0)

        self.filter_spin = QSpinBox()
        self.filter_spin.setRange(0, 5)
        self.filter_spin.setSpecialValueText(self.tr('off'))

        self.image = image
        self.viewer = ImageViewer(self.image, self.image)
        self.low = self.high = None
        self.preprocess()

        self.chan_combo.currentIndexChanged.connect(self.preprocess)
        self.min_combo.currentIndexChanged.connect(self.process)
        self.max_combo.currentIndexChanged.connect(self.process)
        self.filter_spin.valueChanged.connect(self.process)

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr('Channel:')))
        top_layout.addWidget(self.chan_combo)
        top_layout.addWidget(QLabel(self.tr('Minimum:')))
        top_layout.addWidget(self.min_combo)
        top_layout.addWidget(QLabel(self.tr('Maximum:')))
        top_layout.addWidget(self.max_combo)
        top_layout.addWidget(QLabel(self.tr('Filter:')))
        top_layout.addWidget(self.filter_spin)
        top_layout.addStretch()
        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    @staticmethod
    def minmax_dev(patch, mask):
        c = patch[1, 1]
        minimum, maximum, _, _ = cv.minMaxLoc(patch, mask)
        if c < minimum:
            return -1
        if c > maximum:
            return +1
        return 0

    @staticmethod
    def blk_filter(img, radius):
        result = np.zeros_like(img, np.float32)
        rows, cols = result.shape
        block = 2*radius + 1
        for i in range(radius, rows, block):
            for j in range(radius, cols, block):
                result[i-radius:i+radius+1, j-radius:j+radius+1] = np.std(img[i-radius:i+radius+1, j-radius:j+radius+1])
        return cv.normalize(result, None, 0, 127, cv.NORM_MINMAX, cv.CV_8UC1)

    def preprocess(self):
        start = time()
        channel = self.chan_combo.currentIndex()
        if channel == 0:
            img = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        elif channel == 4:
            b, g, r = cv.split(self.image.astype(np.float64))
            img = cv.sqrt(cv.pow(b, 2) + cv.pow(g, 2) + cv.pow(r, 2))
        else:
            img = self.image[:, :, 3 - channel]
        kernel = 3
        border = kernel // 2
        shape = (img.shape[0] - kernel + 1, img.shape[1] - kernel + 1, kernel, kernel)
        strides = 2 * img.strides
        patches = np.lib.stride_tricks.as_strided(img, shape=shape, strides=strides)
        patches = patches.reshape((-1, kernel, kernel))
        mask = np.full((kernel, kernel), 255, dtype=np.uint8)
        mask[border, border] = 0
        output = np.array([self.minmax_dev(patch, mask) for patch in patches]).reshape(shape[:-2])
        output = cv.copyMakeBorder(output, border, border, border, border, cv.BORDER_CONSTANT)
        self.low = output == -1
        self.high = output == +1
        self.process()
        self.info_message.emit(self.tr('Min/Max Deviation = {}'.format(elapsed_time(start))))

    def process(self):
        minmax = np.zeros_like(self.image)
        minimum = self.min_combo.currentIndex()
        maximum = self.max_combo.currentIndex()
        radius = self.filter_spin.value()
        if radius > 0:
            start = time()
            radius += 2
            low = self.blk_filter(self.low, radius)
            high = self.blk_filter(self.high, radius)
            if minimum <= 2:
                minmax[:, :, 2 - minimum] = low
            else:
                minmax = low
            if maximum <= 2:
                minmax[:, :, 2 - maximum] += high
            else:
                minmax += high
            minmax = normalize_mat(minmax)
            self.info_message.emit(self.tr('Min/Max Filter = {}'.format(elapsed_time(start))))
        else:
            if minimum == 0:
                minmax[self.low] = [0, 0, 255]
            elif minimum == 1:
                minmax[self.low] = [0, 255, 0]
            elif minimum == 2:
                minmax[self.low] = [255, 0, 0]
            elif minimum == 3:
                minmax[self.low] = [255, 255, 255]
            if maximum == 0:
                minmax[self.high] = [0, 0, 255]
            elif maximum == 1:
                minmax[self.high] = [0, 255, 0]
            elif maximum == 2:
                minmax[self.high] = [255, 0, 0]
            elif maximum == 3:
                minmax[self.high] = [255, 255, 255]
        self.viewer.update_processed(minmax)
