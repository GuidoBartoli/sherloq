from time import time

import cv2 as cv
import numpy as np
from PySide2.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QSpinBox,
    QLabel)

from tools import ToolWidget
from utility import elapsed_time, create_lut
from viewer import ImageViewer


class EchoWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(EchoWidget, self).__init__(parent)

        params_layout = QHBoxLayout()
        params_layout.addWidget(QLabel(self.tr('Radius:')))
        self.radius_spin = QSpinBox()
        self.radius_spin.setRange(1, 50)
        self.radius_spin.setSuffix(self.tr(' px'))
        self.radius_spin.setValue(2)
        self.radius_spin.valueChanged.connect(self.process)
        params_layout.addWidget(self.radius_spin)

        params_layout.addWidget(QLabel(self.tr('Contrast:')))
        self.contrast_spin = QSpinBox()
        self.contrast_spin.setRange(0, 100)
        self.contrast_spin.setSuffix(self.tr(' %'))
        self.contrast_spin.setValue(85)
        self.contrast_spin.valueChanged.connect(self.process)
        params_layout.addWidget(self.contrast_spin)
        params_layout.addStretch()

        self.image = image
        self.viewer = ImageViewer(self.image, self.image, None)
        self.process()

        main_layout = QVBoxLayout()
        main_layout.addLayout(params_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def process(self):
        start = time()
        kernel = 2*self.radius_spin.value() + 1
        contrast = int(self.contrast_spin.value() / 100 * 255)
        lut = create_lut(0, contrast)
        laplace = []
        for channel in cv.split(self.image):
            deriv = np.fabs(cv.Laplacian(channel, cv.CV_64F, None, kernel))
            deriv = cv.normalize(deriv, None, 0, 255, cv.NORM_MINMAX, cv.CV_8UC1)
            laplace.append(cv.LUT(deriv, lut))
        self.viewer.update_processed(cv.merge(laplace))
        self.info_message.emit('Echo Edge Filter = {}'.format(elapsed_time(start)))
