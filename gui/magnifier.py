import cv2 as cv
import numpy as np
from PySide2.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QRadioButton,
    QSpinBox,
    QCheckBox)

from tools import ToolWidget
from utility import auto_lut, equalize_image
from viewer import ImageViewer


class MagnifierWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(MagnifierWidget, self).__init__(parent)

        self.equalize_radio = QRadioButton(self.tr('Equalization'))
        self.contrast_radio = QRadioButton(self.tr('Auto Contrast'))
        self.retinex_radio = QRadioButton(self.tr('Human Retina'))
        self.centile_spin = QSpinBox()
        self.centile_spin.setRange(0, 100)
        self.centile_spin.setValue(20)
        self.centile_spin.setSuffix(self.tr(' %'))
        self.channel_check = QCheckBox(self.tr('By channel'))
        self.equalize_radio.setChecked(True)
        self.last_radio = self.equalize_radio

        self.image = image
        self.viewer = ImageViewer(self.image, self.image)
        self.change()

        self.viewer.viewChanged.connect(self.process)
        self.equalize_radio.clicked.connect(self.change)
        self.contrast_radio.clicked.connect(self.change)
        self.centile_spin.valueChanged.connect(self.change)
        self.channel_check.stateChanged.connect(self.change)
        self.retinex_radio.clicked.connect(self.change)

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr('Enhancement:')))
        top_layout.addWidget(self.equalize_radio)
        top_layout.addWidget(self.contrast_radio)
        top_layout.addWidget(self.centile_spin)
        top_layout.addWidget(self.channel_check)
        top_layout.addWidget(self.retinex_radio)
        top_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def process(self, rect):
        y1 = rect.top()
        y2 = rect.bottom()
        x1 = rect.left()
        x2 = rect.right()
        roi = self.image[y1:y2, x1:x2]
        if self.equalize_radio.isChecked():
            roi = equalize_image(roi)
            self.last_radio = self.equalize_radio
        elif self.contrast_radio.isChecked():
            centile = self.centile_spin.value() / 200
            if self.channel_check.isChecked():
                roi = cv.merge([cv.LUT(c, auto_lut(c, centile)) for c in cv.split(roi)])
            else:
                roi = cv.LUT(roi, auto_lut(cv.cvtColor(roi, cv.COLOR_BGR2GRAY), centile))
            self.last_radio = self.contrast_radio
        elif self.retinex_radio.isChecked():
            self.last_radio = self.retinex_radio
        else:
            self.last_radio.setChecked(True)
            return
        processed = np.copy(self.image)
        processed[y1:y2, x1:x2] = roi
        self.viewer.update_processed(processed)

    def change(self):
        self.process(self.viewer.get_rect())
