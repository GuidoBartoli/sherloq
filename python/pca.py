import cv2 as cv
import numpy as np
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QPushButton,
    QVBoxLayout,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QProgressDialog,
    QRadioButton)

from tools import ToolWidget
from utility import normalize_mat
from viewer import ImageViewer


class PcaWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(PcaWidget, self).__init__(parent)

        params_layout = QHBoxLayout()
        params_layout.addWidget(QLabel(self.tr('Component:')))
        self.comp_combo = QComboBox()
        self.comp_combo.addItems([self.tr('#{}'.format(i + 1)) for i in range(3)])
        self.comp_combo.setEnabled(False)
        self.comp_combo.currentIndexChanged.connect(self.show)
        params_layout.addWidget(self.comp_combo)

        self.distvect_radio = QRadioButton(self.tr('Vector Distance'))
        self.distvect_radio.setChecked(True)
        self.distvect_radio.setEnabled(False)
        self.distvect_radio.toggled.connect(self.show)
        params_layout.addWidget(self.distvect_radio)
        self.disteq_radio = QRadioButton(self.tr('Equalized Distance'))
        self.disteq_radio.setEnabled(False)
        self.disteq_radio.toggled.connect(self.show)
        params_layout.addWidget(self.disteq_radio)
        self.cross_radio = QRadioButton(self.tr('Cross-Correlation'))
        self.cross_radio.setEnabled(False)
        self.cross_radio.toggled.connect(self.show)
        params_layout.addWidget(self.cross_radio)

        self.process_button = QPushButton(self.tr('Process'))
        self.process_button.setToolTip(self.tr('Start PCA processing (WARNING: this can take a while!)'))
        self.process_button.clicked.connect(self.process)
        params_layout.addWidget(self.process_button)
        params_layout.addStretch()

        self.image = image
        self.components = []
        self.viewer = ImageViewer(self.image, self.image, None)
        main_layout = QVBoxLayout()
        main_layout.addLayout(params_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def process(self):
        self.components = []
        rows, cols, dims = self.image.shape
        bgr = np.reshape(self.image, (rows * cols, dims)).astype(np.float64)
        q = np.empty((0))
        q, x, y = cv.PCACompute2(bgr, q)
        progress = QProgressDialog(self.tr('RGB PCA Projection...'), None, 0, self.image.size, self)
        progress.setWindowModality(Qt.WindowModal)
        counter = 0
        for d in range(dims):
            v = x[d]
            r = q + v * y[d]
            distance = np.zeros((rows, cols), dtype=np.float64)
            cross = np.zeros((rows, cols, dims), dtype=np.float64)
            # TODO: Provare a vettorizzare le operazioni per maggiore velocità
            for i in range(rows):
                for j in range(cols):
                    p = self.image[i, j]
                    distance[i, j] = cv.norm(np.cross(p - q, p - r)) / cv.norm(r - q)
                    cross[i, j] = np.cross(p, v)
                    counter += 1
                    progress.setValue(counter)
            distance_eq = cv.cvtColor(cv.equalizeHist(normalize_mat(distance)), cv.COLOR_GRAY2BGR)
            distance = normalize_mat(distance, to_bgr=True)
            cross = cv.merge([normalize_mat(c) for c in cv.split(cross)]).astype(np.uint8)
            self.components.extend([distance, distance_eq, cross])
        progress.setValue(self.image.size)
        self.comp_combo.setEnabled(True)
        self.distvect_radio.setEnabled(True)
        self.disteq_radio.setEnabled(True)
        self.cross_radio.setEnabled(True)
        self.process_button.setEnabled(False)
        self.show()

    def show(self):
        index = self.comp_combo.currentIndex()
        if self.distvect_radio.isChecked():
            self.viewer.update_processed(self.components[3*index])
        elif self.disteq_radio.isChecked():
            self.viewer.update_processed(self.components[3*index + 1])
        elif self.cross_radio.isChecked():
            self.viewer.update_processed(self.components[3*index + 2])
