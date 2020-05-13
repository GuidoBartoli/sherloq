import cv2 as cv
import numpy as np
from PySide2.QtWidgets import (
    QSpinBox,
    QCheckBox,
    QHBoxLayout,
    QVBoxLayout,
    QLabel)

from utility import create_lut
from viewer import ImageViewer
from widget import ToolWidget


class GradientWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(GradientWidget, self).__init__(parent)
        self.intensity_spin = QSpinBox()
        self.intensity_spin.setRange(-1, 16)
        self.intensity_spin.setSuffix(self.tr(' levels'))
        self.intensity_spin.setSpecialValueText(self.tr('Equalized'))
        self.intensity_spin.setValue(-1)
        self.invert_check = QCheckBox(self.tr('Invert'))
        self.abs_check = QCheckBox(self.tr('Absolute'))

        horiz_layout = QHBoxLayout()
        horiz_layout.addWidget(QLabel(self.tr('Intensity')))
        horiz_layout.addWidget(self.intensity_spin)
        horiz_layout.addWidget(self.invert_check)
        horiz_layout.addWidget(self.abs_check)
        horiz_layout.addStretch()

        self.grad_viewer = ImageViewer(image, image)
        self.image = image
        self.process()
        vert_layout = QVBoxLayout()
        vert_layout.addLayout(horiz_layout)
        vert_layout.addWidget(self.grad_viewer)
        self.setLayout(vert_layout)

        self.intensity_spin.valueChanged.connect(self.process)
        self.invert_check.toggled.connect(self.process)
        self.abs_check.toggled.connect(self.process)

    def process(self):
        intensity = self.intensity_spin.value()
        invert = self.invert_check.isChecked()
        absolute = self.abs_check.isChecked()
        self.intensity_spin.setEnabled(not absolute)
        self.invert_check.setEnabled(not absolute)

        gray = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        diff_x = cv.Scharr(gray, cv.CV_32F, 1, 0)
        diff_y = cv.Scharr(gray, cv.CV_32F, 0, 1)
        diff_z = cv.normalize(np.abs(diff_x + diff_y), None, 0, 255, cv.NORM_MINMAX)
        diff_z = cv.LUT(diff_z.astype(np.uint8), create_lut(0, 90))
        grad_x = +diff_x if invert else -diff_x
        grad_y = +diff_y if invert else -diff_y
        grad_z = diff_z

        if not absolute:
            min_x, max_x, _, _ = cv.minMaxLoc(grad_x)
            lim_x = max(abs(min_x), abs(max_x))
            grad_x = (grad_x / lim_x + 1)*127
            min_y, max_y, _, _ = cv.minMaxLoc(grad_y)
            lim_y = max(abs(min_y), abs(max_y))
            grad_y = (grad_y / lim_y + 1)*127
        else:
            grad_x = cv.normalize(np.abs(grad_x), None, 0, 255, cv.NORM_MINMAX)
            grad_y = cv.normalize(np.abs(grad_y), None, 0, 255, cv.NORM_MINMAX)
        grad_x = grad_x.astype(np.uint8)
        grad_y = grad_y.astype(np.uint8)

        if intensity >= 0 and not absolute:
            max_intensity = self.intensity_spin.maximum()
            low = 127 - (max_intensity - intensity)
            high = 127 + (max_intensity - intensity)
            lut_xy = create_lut(low, high)
            grad_x = cv.LUT(grad_x, lut_xy)
            grad_y = cv.LUT(grad_y, lut_xy)
        else:
            grad_x = cv.equalizeHist(grad_x)
            grad_y = cv.equalizeHist(grad_y)
        gradient = np.stack((grad_z, grad_x, grad_y), axis=2)
        self.grad_viewer.update_processed(gradient)
