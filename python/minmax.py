import cv2 as cv
import numpy as np
from PySide2.QtCore import Qt
from PySide2.QtWidgets import QVBoxLayout, QProgressDialog

from viewer import ImageViewer
from widget import ToolWidget


class MinMaxWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(MinMaxWidget, self).__init__(parent)
        min_max = np.zeros_like(image)
        count = 0
        progress = QProgressDialog(self.tr('Computing min/max deviation...'), None, 0, min_max.size // 3, parent)
        progress.setWindowModality(Qt.WindowModal)
        for i in range(1, min_max.shape[0] - 1):
            for j in range(1, min_max.shape[1] - 1):
                norm0 = cv.norm(image[i, j, :])
                norm1 = [cv.norm(image[i, j - 1, :]), cv.norm(image[i - 1, j, :]),
                         cv.norm(image[i, j + 1, :]), cv.norm(image[i + 1, j, :])]
                max_norm = max(norm1)
                min_norm = min(norm1)
                if norm0 > max_norm:
                    min_max[i, j, :] = [0, 0, 255]
                elif norm0 < min_norm:
                    min_max[i, j, :] = [0, 255, 0]
                progress.setValue(count)
                count += 1
        viewer = ImageViewer(image, min_max)
        layout = QVBoxLayout()
        layout.addWidget(viewer)
        self.setLayout(layout)
