from time import time

import cv2 as cv
import numpy as np
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QComboBox,
    QCheckBox,
    QLabel,
    QProgressDialog)

from tools import ToolWidget
from utility import elapsed_time
from viewer import ImageViewer


class MinMaxWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(MinMaxWidget, self).__init__(parent)

        params_layout = QHBoxLayout()
        params_layout.addWidget(QLabel(self.tr('Channel:')))
        self.chan_combo = QComboBox()
        self.chan_combo.addItems(
            [self.tr('Luminance'), self.tr('Red'), self.tr('Green'), self.tr('Blue'), self.tr('L2-Norm')])
        params_layout.addWidget(self.chan_combo)

        colors = [self.tr('Red'), self.tr('Green'), self.tr('Blue'), self.tr('White'), self.tr('Black')]
        params_layout.addWidget(QLabel(self.tr('Minimum:')))
        self.min_combo = QComboBox()
        self.min_combo.addItems(colors)
        self.min_combo.setCurrentIndex(1)
        params_layout.addWidget(self.min_combo)

        params_layout.addWidget(QLabel(self.tr('Maximum:')))
        self.max_combo = QComboBox()
        self.max_combo.addItems(colors)
        self.max_combo.setCurrentIndex(0)
        params_layout.addWidget(self.max_combo)

        self.accurate_check = QCheckBox(self.tr('Accurate'))
        params_layout.addWidget(self.accurate_check)
        params_layout.addStretch()

        self.image = image
        self.viewer = ImageViewer(self.image, self.image)
        self.process()

        main_layout = QVBoxLayout()
        main_layout.addLayout(params_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

        self.chan_combo.currentIndexChanged.connect(self.process)
        self.min_combo.currentIndexChanged.connect(self.process)
        self.max_combo.currentIndexChanged.connect(self.process)
        self.accurate_check.stateChanged.connect(self.process)

    def process(self):
        start = time()
        channel = self.chan_combo.currentIndex()
        if channel == 0:
            original = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        elif channel == 4:
            b, g, r = cv.split(self.image.astype(np.float64))
            original = cv.sqrt(cv.pow(b, 2) + cv.pow(g, 2) + cv.pow(r, 2))
        else:
            original = self.image[:, :, 3 - channel]
        accurate = self.accurate_check.isChecked()
        if not accurate:
            kernel = cv.getStructuringElement(cv.MORPH_RECT, (3, 3))
            low = original == cv.erode(original, kernel)
            high = original == cv.dilate(original, kernel)
        else:
            progress = QProgressDialog(self.tr('Min/Max Deviation (accurate)...'), None, 0, original.size, self)
            progress.setWindowModality(Qt.WindowModal)
            low = np.full_like(original, False, dtype=bool)
            high = np.full_like(original, False, dtype=bool)
            mask = np.ones((3, 3), dtype=np.uint8)
            mask[1, 1] = 0
            rows, cols = original.shape
            counter = 0
            for i in range(1, rows - 1):
                for j in range(1, cols - 1):
                    progress.setValue(counter)
                    window = original[i-1:i+2, j-1:j+2]
                    minval, maxval, _, _ = cv.minMaxLoc(window, mask)
                    if original[i, j] < minval:
                        low[i, j] = True
                    elif original[i, j] > maxval:
                        high[i, j] = True
                    counter += 1
            progress.setValue(original.size)
        minmax = np.zeros_like(self.image)
        minimum = self.min_combo.currentIndex()
        maximum = self.max_combo.currentIndex()

        if minimum == 0:
            minmax[low] = [0, 0, 255]
        elif minimum == 1:
            minmax[low] = [0, 255, 0]
        elif minimum == 2:
            minmax[low] = [255, 0, 0]
        elif minimum == 3:
            minmax[low] = [255, 255, 255]

        if maximum == 0:
            minmax[high] = [0, 0, 255]
        elif maximum == 1:
            minmax[high] = [0, 255, 0]
        elif maximum == 2:
            minmax[high] = [255, 0, 0]
        elif maximum == 3:
            minmax[high] = [255, 255, 255]

        self.viewer.update_processed(minmax)
        self.info_message.emit(self.tr('Min/Max Deviation = {}'.format(elapsed_time(start))))
