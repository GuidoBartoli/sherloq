import cv2 as cv
import numpy as np
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QSpinBox,
    QComboBox,
    QCheckBox,
    QHBoxLayout,
    QVBoxLayout,
    QProgressDialog,
    QPushButton,
    QLabel)

from tools import ToolWidget
from utility import create_lut, norm_mat, equalize_img, elapsed_time, compute_hist
from viewer import ImageViewer


class ContrastWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(ContrastWidget, self).__init__(parent)

        rows0, cols0, _ = image.shape
        block = 128 if min(rows0, cols0) / 128 > 10 else 64
        color = cv.copyMakeBorder(image, 0, rows0 % block, 0, cols0 % block, cv.BORDER_CONSTANT)
        gray = cv.cvtColor(color, cv.COLOR_BGR2GRAY)
        rows, cols = gray.shape

        kx, ky = cv.getDerivKernels(1, 1, 1)
        bd, gd, rd = [cv.sepFilter2D(c, cv.CV_32F, kx, ky) for c in cv.split(color)]
        tri = (np.abs(gd - rd) + np.abs(gd - bd) + np.abs(rd - bd)) / 3
        avg = (np.abs(bd) + np.abs(gd) + np.abs(rd)) / 3

        window = np.arange(256).astype(np.float32)
        cutoff = 8
        window[:cutoff] = (1 - np.cos(np.pi * window[:cutoff] / cutoff)) / 2
        window[-cutoff:] = (1 + np.cos(np.pi * (window[-cutoff:] + cutoff - 255) / cutoff)) / 2
        window[cutoff:-cutoff] = 1
        weight = ((np.arange(256) - 128) / 128) ** 2

        output = np.zeros(((rows // block) + 1, (cols // block) + 1), np.float32)
        progress = QProgressDialog(self.tr('Detecting enhancements...'), None, 0, output.shape[0] - 1, parent)
        progress.setWindowModality(Qt.WindowModal)
        max_err = 0.185
        max_sim = 0.75
        p = 0
        for i in range(0, rows, block):
            for j in range(0, cols, block):
                hist = compute_hist(gray[i:i+block, j:j+block]) * window
                hist = cv.normalize(hist, None, 0, 1, cv.NORM_MINMAX)
                dft = np.fft.fftshift(cv.dft(hist, flags=cv.DFT_COMPLEX_OUTPUT))
                mag = cv.magnitude(dft[:, :, 0], dft[:, :, 1])
                mag = cv.normalize(mag, None, 0, 1, cv.NORM_MINMAX).flatten()

                diff = 0
                for k in range(2, 254):
                    yl = 2 * hist[k - 1] - hist[k - 2]
                    yr = 2 * hist[k + 1] - hist[k + 2]
                    d = abs(hist[k] - (yl + yr) / 2)
                    if d > diff:
                        diff = d

                ed = np.sum(mag)
                if ed == 0:
                    error = 0
                else:
                    en = np.sum(mag * weight)
                    error = en / ed
                    if error > max_err:
                        error = 1
                    else:
                        error /= max_err
                    error *= np.sqrt(diff)

                avg_m = np.mean(avg[i:i+block, j:j+block])
                if avg_m == 0:
                    chsim = 0
                else:
                    tri_m = np.mean(tri[i:i+block, j:j+block])
                    chsim = tri_m / avg_m
                    if chsim > max_sim:
                        chsim = 1
                    else:
                        chsim /= max_sim
                output[i // block, j // block] = error * chsim
            progress.setValue(p)
            p += 1

        output = cv.medianBlur(cv.convertScaleAbs(output, None, 255), 3)
        output = cv.resize(output, None, None, block, block, cv.INTER_NEAREST)[:rows0, :cols0]
        viewer = ImageViewer(image, cv.cvtColor(output, cv.COLOR_GRAY2BGR))
        layout = QVBoxLayout()
        layout.addWidget(viewer)
        self.setLayout(layout)

    def process(self):
        pass
