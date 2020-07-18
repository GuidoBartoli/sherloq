import os
from time import time

import cv2 as cv
import numpy as np
from joblib import load
import sewar
from PySide2.QtCore import QCoreApplication, Qt
from PySide2.QtWidgets import (
    QComboBox,
    QLabel,
    QHBoxLayout,
    QPushButton,
    QVBoxLayout,
    QProgressDialog,
    QCheckBox,
    QGridLayout,
    QMessageBox)

from tools import ToolWidget
from utility import modify_font, norm_mat, gray_to_bgr
from viewer import ImageViewer


def ssim(a, b, maximum=255):
    c1 = (0.01 * maximum) ** 2
    c2 = (0.03 * maximum) ** 2
    k = (11, 11)
    s = 1.5
    a2 = a ** 2
    b2 = b ** 2
    ab = a * b
    mu_a = cv.GaussianBlur(a, k, s)
    mu_b = cv.GaussianBlur(b, k, s)
    mu_a2 = mu_a ** 2
    mu_b2 = mu_b ** 2
    mu_ab = mu_a * mu_b
    s_a2 = cv.GaussianBlur(a2, k, s) - mu_a2
    s_b2 = cv.GaussianBlur(b2, k, s) - mu_b2
    s_ab = cv.GaussianBlur(ab, k, s) - mu_ab
    t1 = 2 * mu_ab + c1
    t2 = 2 * s_ab + c2
    t3 = t1 * t2
    t1 = mu_a2 + mu_b2 + c1
    t2 = s_a2 + s_b2 + c2
    t1 *= t2
    s_map = cv.divide(t3, t1)
    return cv.mean(s_map)[0]


def get_metrics(pristine, distorted, normalized=False, extended=False):
    # Matrix precomputation
    x0 = pristine.astype(np.float64)
    y0 = distorted.astype(np.float64)
    if normalized:
        x0 /= 255
        y0 /= 255
    x2 = np.sum(np.square(x0))
    y2 = np.sum(np.square(y0))
    xs = np.sum(x0)
    e = x0 - y0
    maximum = 255 if not normalized else 1
    # Feature vector initialization
    m = np.zeros(8 if not extended else 16)
    # Mean Square Error (MSE)
    m[0] = np.mean(np.square(e))
    # Peak to Signal Noise Ratio (PSNR)
    m[1] = 20 * np.log10(maximum / np.sqrt(m[0])) if m[0] > 0 else -1
    # Normalized Cross-Correlation (NCC)
    m[2] = np.sum(x0 * y0) / x2 if x2 > 0 else -1
    # Average Difference (AD)
    m[3] = np.mean(e)
    # Structural Content (SC)
    m[4] = x2 / y2 if y2 > 0 else -1
    # Maximum Difference (MD)
    m[5] = np.max(e)
    # Normalized Absolute Error (NAE)
    m[6] = np.sum(np.abs(e)) / xs if xs > 0 else -1
    # Structural Similarity (SSIM)
    m[7] = ssim(x0, y0, maximum)
    if extended:
        pass
    return m


def get_features(image, windows=4, levels=4, normalized=True, extended=False):
    metrics = 8 if not extended else 16
    f = np.zeros(windows * levels * metrics)
    index = 0
    for w in range(windows):
        k = 2 * (w + 1) + 1
        previous = image
        for _ in range(levels):
            filtered = cv.medianBlur(previous, k)
            f[index:index+metrics] = get_metrics(previous, filtered, normalized, extended)
            index += metrics
            previous = filtered
    return f


class MedianWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(MedianWidget, self).__init__(parent)

        self.block_combo = QComboBox()
        self.block_combo.addItems(['16', '32', '64', '128', '256'])
        self.block_combo.setCurrentIndex(2)
        self.block_combo.setToolTip(self.tr('Size of analyzed blocks'))
        self.multi_check = QCheckBox(self.tr('Multi-scale'))
        self.multi_check.setChecked(True)
        self.process_button = QPushButton(self.tr('Process'))
        self.prob_label = QLabel()

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr('Block size:')))
        top_layout.addWidget(self.block_combo)
        # top_layout.addWidget(self.multi_check)
        top_layout.addWidget(self.process_button)
        top_layout.addWidget(self.prob_label)
        top_layout.addStretch()

        self.image = image
        self.gray = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        self.viewer = ImageViewer(self.image, self.image)
        self.canceled = False

        self.process_button.clicked.connect(self.process)
        self.block_combo.currentIndexChanged.connect(self.reset)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def reset(self):
        self.process_button.setEnabled(True)

    def cancel(self):
        self.canceled = True

    def process(self):
        rows0, cols0 = self.gray.shape
        block = int(self.block_combo.currentText())
        padded = cv.copyMakeBorder(self.gray, 0, block - rows0 % block, 0, block - cols0 % block, cv.BORDER_CONSTANT)
        rows, cols = padded.shape
        prob = np.zeros(((rows // block) + 1, (cols // block) + 1))

        # if self.multi_check.isChecked():
        #     model = load('models/median_multi.mdl')
        # else:
        #     model = load('models/median_single.mdl')
        model = load('models/fpmw_block.mdl')
        columns = model._features_count
        if columns == 8:
            levels = 1
            windows = 1
        elif columns == 24:
            levels = 3
            windows = 1
        elif columns == 96:
            levels = 3
            windows = 4
        elif columns == 128:
            levels = 4
            windows = 4
        else:
            return
        limit = model.best_ntree_limit if hasattr(model, 'best_ntree_limit') else None

        if not self.prob_label.text():
            self.prob_label.setText(self.tr('Processing, please wait...'))
            modify_font(self.prob_label, italic=True)
            QCoreApplication.processEvents()
            features = np.reshape(get_features(self.gray, levels, windows), (1, columns))
            p = model.predict_proba(features, ntree_limit=limit)[0, 1]
            self.prob_label.setText(self.tr('Filtering probability (whole image) = {:.2f}%'.format(p * 100)))
            modify_font(self.prob_label, italic=False, bold=True)
            QCoreApplication.processEvents()

        progress = QProgressDialog(self.tr('Detecting median filter...'), self.tr('Cancel'), 0, prob.size, self)
        progress.canceled.connect(self.cancel)
        progress.setWindowModality(Qt.WindowModal)
        p = 0
        for i in range(0, rows, block):
            for j in range(0, cols, block):
                roi = padded[i:i+block, j:j+block]
                if np.var(roi) < 5:
                    prob[i // block, j // block] = 0
                else:
                    x = np.reshape(get_features(roi, levels, windows), (1, columns))
                    prob[i // block, j // block] = model.predict_proba(x, ntree_limit=limit)[0, 1]
                if self.canceled:
                    self.canceled = False
                    progress.close()
                    return
                progress.setValue(p)
                p += 1
        progress.setValue(prob.size)

        prob = cv.convertScaleAbs(prob, None, 255)
        prob = cv.resize(prob, None, None, block, block, cv.INTER_NEAREST)
        prob = gray_to_bgr(prob[:rows0, :cols0])
        self.viewer.update_processed(prob)
        self.process_button.setEnabled(False)
