import cv2 as cv
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QSpinBox,
    QSlider,
    QLabel)

from utility import compress_jpeg
from viewer import ImageViewer
from widget import ToolWidget


class ElaWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(ElaWidget, self).__init__(parent)

        self.quality_slider = QSlider(Qt.Horizontal)
        self.quality_slider.setRange(0, 100)
        self.quality_slider.setSingleStep(1)
        self.quality_slider.setPageStep(5)
        self.quality_label = QLabel()
        self.scale_spin = QSpinBox()
        self.scale_spin.setRange(0, 100)
        self.default_button = QPushButton(self.tr('Default'))

        params_layout = QHBoxLayout()
        params_layout.addWidget(QLabel(self.tr("Quality:")))
        params_layout.addWidget(self.quality_slider)
        params_layout.addWidget(self.quality_label)
        params_layout.addWidget(QLabel(self.tr("Scale:")))
        params_layout.addWidget(self.scale_spin)
        params_layout.addStretch()
        params_layout.addWidget(self.default_button)

        self.ela_viewer = ImageViewer(image, image, None)
        self.image = image
        self.default()
        self.process()

        layout = QVBoxLayout()
        layout.addLayout(params_layout)
        layout.addWidget(self.ela_viewer)
        self.setLayout(layout)
        self.quality_slider.valueChanged.connect(self.process)
        self.scale_spin.valueChanged.connect(self.process)
        self.default_button.clicked.connect(self.default)

    def process(self):
        quality = self.quality_slider.value()
        scale = self.scale_spin.value()
        compressed = compress_jpeg(self.image, quality)
        ela = cv.convertScaleAbs(cv.subtract(compressed, self.image), None, scale)
        self.ela_viewer.update_processed(ela)
        self.quality_label.setText('{}%'.format(self.quality_slider.value()))

    def default(self):
        self.quality_slider.setValue(75)
        self.scale_spin.setValue(20)
