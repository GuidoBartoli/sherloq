import math

import cv2 as cv
import numpy as np
from PySide2.QtWidgets import (
    QTableWidgetItem,
    QTableWidget,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QCheckBox,
    QLabel,
    QRadioButton)

from tools import ToolWidget
from utility import normalize_mat, modify_font, load_image, desaturate
from viewer import ImageViewer


class ComparisonWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(ComparisonWidget, self).__init__(parent)

        load_button = QPushButton(self.tr('Load reference...'))
        self.file_label = QLabel()
        modify_font(self.file_label, bold=True)
        self.comp_label = QLabel(self.tr('Comparison:'))
        self.normal_radio = QRadioButton(self.tr('Normal'))
        self.normal_radio.setChecked(True)
        self.diff_radio = QRadioButton(self.tr('Difference'))
        self.gray_check = QCheckBox(self.tr('Grayscale'))
        self.equal_check = QCheckBox(self.tr('Equalized'))
        self.map_radio = QRadioButton(self.tr('SSIM Map'))
        self.last_radio = self.normal_radio

        self.comp_label.setEnabled(False)
        self.normal_radio.setEnabled(False)
        self.diff_radio.setEnabled(False)
        self.gray_check.setEnabled(False)
        self.equal_check.setEnabled(False)
        self.map_radio.setEnabled(False)

        self.evidence = image
        self.reference = np.full_like(image, 127)
        self.difference = np.full_like(image, 127)
        self.equalized = np.full_like(image, 127)
        self.indexmap = np.full_like(image, 127)
        self.evidence_viewer = ImageViewer(self.evidence, None, self.tr('Evidence'))
        self.reference_viewer = ImageViewer(self.reference, None, self.tr('Reference'))

        self.table_widget = QTableWidget(5, 2)
        self.table_widget.setHorizontalHeaderLabels([self.tr('Index'), self.tr('Value')])
        self.table_widget.setItem(0, 0, QTableWidgetItem(self.tr('MSE')))
        self.table_widget.setItem(1, 0, QTableWidgetItem(self.tr('COVAR')))
        self.table_widget.setItem(2, 0, QTableWidgetItem(self.tr('PSNR')))
        self.table_widget.setItem(3, 0, QTableWidgetItem(self.tr('SSIM')))
        self.table_widget.setItem(4, 0, QTableWidgetItem(self.tr('CORR')))
        for i in range(self.table_widget.rowCount()):
            modify_font(self.table_widget.item(i, 0), bold=True)
        self.table_widget.setEnabled(False)

        load_button.clicked.connect(self.load)
        self.normal_radio.toggled.connect(self.change)
        self.diff_radio.toggled.connect(self.change)
        self.gray_check.stateChanged.connect(self.change)
        self.equal_check.stateChanged.connect(self.change)
        self.map_radio.toggled.connect(self.change)
        self.evidence_viewer.view_changed.connect(self.reference_viewer.change_view)
        self.reference_viewer.view_changed.connect(self.evidence_viewer.change_view)

        top_layout = QHBoxLayout()
        top_layout.addWidget(load_button)
        top_layout.addWidget(self.file_label)
        top_layout.addStretch()
        top_layout.addWidget(self.comp_label)
        top_layout.addWidget(self.normal_radio)
        top_layout.addWidget(self.map_radio)
        top_layout.addWidget(self.diff_radio)
        top_layout.addWidget(self.gray_check)
        top_layout.addWidget(self.equal_check)

        index_layout = QVBoxLayout()
        index_label = QLabel(self.tr('Image Quality Assessment'))
        modify_font(index_label, bold=True)
        index_layout.addWidget(index_label)
        index_layout.addWidget(self.table_widget)

        center_layout = QHBoxLayout()
        center_layout.addWidget(self.evidence_viewer)
        center_layout.addWidget(self.reference_viewer)
        center_layout.addLayout(index_layout)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addLayout(center_layout)
        self.setLayout(main_layout)

    def load(self):
        filename, basename, image = load_image(self)
        if filename is None:
            return
        if image.shape != self.evidence.shape:
            QMessageBox.critical(self, self.tr('Error'), self.tr('Evidence and reference must have the same size!'))
            return
        self.file_label.setText(basename)
        self.reference = image
        self.difference = normalize_mat(cv.absdiff(self.evidence, self.reference))
        self.equalized = cv.merge([cv.equalizeHist(c) for c in cv.split(self.difference)])

        x = cv.cvtColor(self.evidence, cv.COLOR_BGR2GRAY).astype(np.float32)
        y = cv.cvtColor(self.reference, cv.COLOR_BGR2GRAY).astype(np.float32)
        mse = self.mse(x, y)
        covar = self.covar(x, y)
        psnr = self.psnr(mse)
        ssim, self.indexmap = self.ssim(x, y)
        corr = self.corr(x, y)
        self.table_widget.setItem(0, 1, QTableWidgetItem('{:.4f}'.format(mse)))
        self.table_widget.setItem(1, 1, QTableWidgetItem('{:.4f}'.format(covar)))
        if psnr > 0:
            self.table_widget.setItem(2, 1, QTableWidgetItem('{:.2f} dB'.format(psnr)))
        else:
            self.table_widget.setItem(2, 1, QTableWidgetItem('+' + u'\u221e' + ' dB'))
        self.table_widget.setItem(3, 1, QTableWidgetItem('{:.4f}'.format(ssim)))
        self.table_widget.setItem(4, 1, QTableWidgetItem('{:.4f}'.format(corr)))
        self.table_widget.setEnabled(True)

        self.comp_label.setEnabled(True)
        self.normal_radio.setEnabled(True)
        self.diff_radio.setEnabled(True)
        self.gray_check.setEnabled(True)
        self.equal_check.setEnabled(True)
        self.map_radio.setEnabled(True)
        self.change()

    def change(self):
        self.gray_check.setEnabled(False)
        self.equal_check.setEnabled(False)
        if self.normal_radio.isChecked():
            self.reference_viewer.update_original(self.reference)
            self.last_radio = self.normal_radio
        elif self.diff_radio.isChecked():
            self.gray_check.setEnabled(True)
            self.equal_check.setEnabled(True)
            if self.equal_check.isChecked():
                result = self.equalized
            else:
                result = self.difference
            if self.gray_check.isChecked():
                result = desaturate(result)
            self.reference_viewer.update_original(result)
            self.last_radio = self.diff_radio
        elif self.map_radio.isChecked():
            self.reference_viewer.update_original(self.indexmap)
            self.last_radio = self.map_radio
        else:
            self.last_radio.setChecked(True)

    @staticmethod
    def ssim(x, y):
        c1 = 6.5025
        c2 = 58.5225
        k = (11, 11)
        s = 1.5
        x2 = x ** 2
        y2 = y ** 2
        xy = x * y
        mu_x = cv.GaussianBlur(x, k, s)
        mu_y = cv.GaussianBlur(y, k, s)
        mu_x2 = mu_x ** 2
        mu_y2 = mu_y ** 2
        mu_xy = mu_x * mu_y
        s_x2 = cv.GaussianBlur(x2, k, s) - mu_x2
        s_y2 = cv.GaussianBlur(y2, k, s) - mu_y2
        s_xy = cv.GaussianBlur(xy, k, s) - mu_xy
        t1 = 2 * mu_xy + c1
        t2 = 2 * s_xy + c2
        t3 = t1 * t2
        t1 = mu_x2 + mu_y2 + c1
        t2 = s_x2 + s_y2 + c2
        t1 *= t2
        indexmap = cv.divide(t3, t1)
        ssim = cv.mean(indexmap)[0]
        return ssim, 255 - normalize_mat(indexmap, to_bgr=True)

    @staticmethod
    def corr(x, y):
        return np.corrcoef(x, y)[0, 1]

    @staticmethod
    def covar(x, y):
        return np.std(cv.absdiff(x, y))

    @staticmethod
    def mse(x, y):
        return cv.mean(cv.pow(x - y, 2))[0]

    @staticmethod
    def psnr(mse):
        k = math.sqrt(mse)
        if k == 0:
            return -1
        return 20 * math.log10((255**2) / k)

