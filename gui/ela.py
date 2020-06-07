from time import time

import cv2 as cv
from PySide2.QtWidgets import (
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QCheckBox,
    QSpinBox,
    QLabel)

from jpeg import compress_jpeg
from tools import ToolWidget
from utility import elapsed_time, equalize_img
from viewer import ImageViewer


class ElaWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(ElaWidget, self).__init__(parent)

        self.quality_spin = QSpinBox()
        self.quality_spin.setRange(0, 100)
        self.quality_spin.setSuffix(self.tr(' %'))

        self.scale_spin = QSpinBox()
        self.scale_spin.setRange(1, 100)

        self.equalize_check = QCheckBox(self.tr('Equalized'))
        self.gray_check = QCheckBox(self.tr('Grayscale'))
        default_button = QPushButton(self.tr('Default'))

        self.quality_spin.valueChanged.connect(self.process)
        self.scale_spin.valueChanged.connect(self.process)
        self.equalize_check.stateChanged.connect(self.process)
        self.gray_check.stateChanged.connect(self.process)
        default_button.clicked.connect(self.default)

        params_layout = QHBoxLayout()
        params_layout.addWidget(QLabel(self.tr('Quality:')))
        params_layout.addWidget(self.quality_spin)
        params_layout.addWidget(QLabel(self.tr('Scale:')))
        params_layout.addWidget(self.scale_spin)
        params_layout.addWidget(self.equalize_check)
        params_layout.addWidget(self.gray_check)
        params_layout.addWidget(default_button)
        params_layout.addStretch()

        self.image = image
        self.viewer = ImageViewer(self.image, self.image)
        self.default()
        self.process()

        main_layout = QVBoxLayout()
        main_layout.addLayout(params_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def process(self):
        start = time()
        quality = self.quality_spin.value()
        scale = self.scale_spin.value()
        equalize = self.equalize_check.isChecked()
        grayscale = self.gray_check.isChecked()
        self.scale_spin.setEnabled(not equalize)
        compressed = compress_jpeg(self.image, quality)
        if not equalize:
            # TODO: Provare a replicare il risultato di FotoForensic dove si vedono di più i blocchi JPEG
            ela = cv.convertScaleAbs(cv.subtract(compressed, self.image), None, scale)
        else:
            ela = equalize_img(cv.absdiff(compressed, self.image))
        if grayscale:
            ela = cv.cvtColor(cv.cvtColor(ela, cv.COLOR_BGR2GRAY), cv.COLOR_GRAY2BGR)
        self.viewer.update_processed(ela)
        self.info_message.emit(self.tr('Error Level Analysis = {}'.format(elapsed_time(start))))

    def default(self):
        self.quality_spin.setValue(75)
        self.scale_spin.setValue(20)
