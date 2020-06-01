import cv2 as cv
import numpy as np
from PySide2.QtWidgets import (
    QHBoxLayout,
    QRadioButton,
    QCheckBox,
    QLabel,
    QVBoxLayout,
    QSpinBox)

from tools import ToolWidget
from utility import normalize_mat
from viewer import ImageViewer


class FrequencyWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(FrequencyWidget, self).__init__(parent)

        self.ampl_radio = QRadioButton(self.tr('Amplitude'))
        self.phase_radio = QRadioButton(self.tr('Phase'))
        self.dct_radio = QRadioButton(self.tr('DCT Map'))
        self.last_radio = self.ampl_radio
        self.thr_spin = QSpinBox()
        self.thr_spin.setRange(0, 255)
        self.thr_spin.setSpecialValueText(self.tr('Off'))
        self.ratio_label = QLabel()
        self.filter_check = QCheckBox(self.tr('Filter'))

        self.ampl_radio.clicked.connect(self.process)
        self.phase_radio.clicked.connect(self.process)
        self.dct_radio.clicked.connect(self.process)
        self.thr_spin.valueChanged.connect(self.process)
        self.filter_check.stateChanged.connect(self.process)

        gray = cv.cvtColor(image, cv.COLOR_BGR2GRAY)
        rows, cols = gray.shape
        height = cv.getOptimalDFTSize(rows)
        width = cv.getOptimalDFTSize(cols)
        padded = cv.copyMakeBorder(gray, 0, height - rows, 0, width - cols, cv.BORDER_CONSTANT).astype(np.float32)
        planes = cv.merge([padded, np.zeros_like(padded)])
        dft = cv.split(np.fft.fftshift(cv.dft(planes)))
        mag, phase = cv.cartToPolar(dft[0], dft[1])
        dct = cv.dct(padded)
        self.result = [normalize_mat(img) for img in [cv.log(mag), phase, cv.log(dct)]]

        self.image = image
        self.viewer = ImageViewer(self.image, None)
        self.process()

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr('Coefficients:')))
        top_layout.addWidget(self.ampl_radio)
        top_layout.addWidget(self.phase_radio)
        top_layout.addWidget(self.dct_radio)
        top_layout.addWidget(self.filter_check)
        top_layout.addStretch()
        top_layout.addWidget(QLabel(self.tr('Threshold:')))
        top_layout.addWidget(self.thr_spin)
        top_layout.addWidget(self.ratio_label)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def process(self):
        if self.ampl_radio.isChecked():
            output = self.result[0]
            self.last_radio = self.ampl_radio
        elif self.phase_radio.isChecked():
            output = self.result[1]
            self.last_radio = self.phase_radio
        elif self.dct_radio.isChecked():
            output = self.result[2]
            self.last_radio = self.dct_radio
        else:
            self.last_radio.setChecked(True)
            return
        if self.filter_check.isChecked():
            output = cv.medianBlur(output, 3)
        thr = self.thr_spin.value()
        if thr > 0:
            _, output = cv.threshold(output, thr, 0, cv.THRESH_TOZERO)
            zeros = (1.0 - cv.countNonZero(output) / output.size) * 100
        else:
            zeros = 0
        self.ratio_label.setText('(masked = {:.1f}%)'.format(zeros))
        self.viewer.update_original(cv.cvtColor(output, cv.COLOR_GRAY2BGR))
