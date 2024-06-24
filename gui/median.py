import cv2 as cv
import numpy as np
import xgboost as xgb
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDoubleSpinBox,
    QMessageBox,
    QSpinBox,
    QLabel,
    QHBoxLayout,
    QPushButton,
    QVBoxLayout,
    QProgressDialog,
    QCheckBox,
)
from joblib import load

from tools import ToolWidget
from utility import modify_font, pad_image
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


def get_metrics(pristine, distorted):
    # Matrix precomputation
    x0 = pristine.astype(np.float64)
    y0 = distorted.astype(np.float64)
    x2 = np.sum(np.square(x0))
    y2 = np.sum(np.square(y0))
    xs = np.sum(x0)
    e = x0 - y0
    maximum = 255
    # Feature vector initialization
    m = np.zeros(8)
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
    return m


def get_features(image, windows, levels):
    metrics = 8
    f = np.zeros(windows * levels * metrics)
    index = 0
    for w in range(windows):
        k = 2 * (w + 1) + 1
        previous = image
        for _ in range(levels):
            filtered = cv.medianBlur(previous, k)
            f[index : index + metrics] = get_metrics(previous, filtered)
            index += metrics
            previous = filtered
    return f


class MedianWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(MedianWidget, self).__init__(parent)

        self.variance_spin = QSpinBox()
        self.variance_spin.setRange(0, 100)
        self.variance_spin.setValue(5)
        self.threshold_spin = QDoubleSpinBox()
        self.threshold_spin.setRange(0, 1)
        self.threshold_spin.setValue(0.4)
        self.threshold_spin.setSingleStep(0.01)
        self.showprob_check = QCheckBox(self.tr("Probability map"))
        self.filter_check = QCheckBox(self.tr("Speckle filter"))
        self.filter_check.setChecked(True)
        self.process_button = QPushButton(self.tr("Process"))
        self.avgprob_label = QLabel(self.tr('Press "Process" to start'))

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr("Min variance:")))
        top_layout.addWidget(self.variance_spin)
        top_layout.addWidget(QLabel(self.tr("Threshold:")))
        top_layout.addWidget(self.threshold_spin)
        top_layout.addWidget(self.showprob_check)
        top_layout.addWidget(self.filter_check)
        top_layout.addWidget(self.process_button)
        # top_layout.addWidget(self.avgprob_label)
        top_layout.addStretch()

        self.image = image
        self.viewer = ImageViewer(self.image, self.image)
        self.gray = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        self.prob = self.var = None
        self.block = 64
        self.canceled = False
        self.modelfile = f"models/median_b{self.block}.json"

        self.process_button.clicked.connect(self.prepare)
        self.variance_spin.valueChanged.connect(self.process)
        self.threshold_spin.valueChanged.connect(self.process)
        self.showprob_check.stateChanged.connect(self.process)
        self.filter_check.stateChanged.connect(self.process)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def prepare(self):
        booster = xgb.Booster()
        try:
            booster.load_model(self.modelfile)
        except xgb.core.XGBoostError:
            QMessageBox.critical(
                self,
                self.tr("Error"),
                self.tr(f'Unable to load model ("{modelfile}")!'),
            )
            return
        columns = booster.num_features()
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
            QMessageBox.critical(
                self, self.tr("Error"), self.tr("Unknown model format!")
            )
            return

        padded = pad_image(self.gray, self.block)
        rows, cols = padded.shape
        self.prob = np.zeros(((rows // self.block) + 1, (cols // self.block) + 1))
        self.var = np.zeros_like(self.prob)
        progress = QProgressDialog(
            self.tr("Detecting median filter..."),
            self.tr("Cancel"),
            0,
            self.prob.size,
            self,
        )
        progress.canceled.connect(self.cancel)
        progress.setWindowModality(Qt.WindowModal)
        k = 0
        self.canceled = False
        for i in range(0, rows, self.block):
            for j in range(0, cols, self.block):
                roi = padded[i : i + self.block, j : j + self.block]
                x = xgb.DMatrix(
                    np.reshape(get_features(roi, levels, windows), (1, columns))
                )
                y = booster.predict(x)[0]
                ib = i // self.block
                jb = j // self.block
                self.var[ib, jb] = np.var(roi)
                self.prob[ib, jb] = y
                if self.canceled:
                    self.prob = self.var = None
                    progress.close()
                    return
                progress.setValue(k)
                k += 1
        progress.close()
        self.process()

    def cancel(self):
        self.canceled = True

    def process(self):
        if self.prob is None:
            return
        mask = self.var < self.variance_spin.value()
        if self.filter_check.isChecked():
            prob = cv.medianBlur(self.prob.astype(np.float32), 3)
        else:
            prob = self.prob.astype(np.float32)
        if self.showprob_check.isChecked():
            output = np.repeat(prob[:, :, np.newaxis], 3, axis=2)
            output[mask] = 0
        else:
            thr = self.threshold_spin.value()
            output = np.zeros((prob.shape[0], prob.shape[1], 3))
            blue, green, red = cv.split(output)
            blue[mask] = 1
            green[prob < thr] = 1
            green[mask] = 0
            red[prob >= thr] = 1
            red[mask] = 0
            output = cv.merge((blue, green, red))
        output = cv.convertScaleAbs(output, None, 255)
        output = cv.resize(output, None, None, self.block, self.block, cv.INTER_LINEAR)
        self.viewer.update_processed(
            np.copy(output[: self.image.shape[0], : self.image.shape[1]])
        )
        avgprob = cv.mean(prob, 1 - mask.astype(np.uint8))[0] * 100
        self.avgprob_label.setText(self.tr(f"Average = {avgprob:.2f}%"))
        modify_font(self.avgprob_label, italic=False, bold=True)
        self.process_button.setEnabled(False)
