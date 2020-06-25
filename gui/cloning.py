import cv2 as cv
import numpy as np
from PySide2.QtWidgets import (
    QPushButton,
    QSpinBox,
    QCheckBox,
    QComboBox,
    QLabel,
    QHBoxLayout,
    QVBoxLayout)

from tools import ToolWidget
from viewer import ImageViewer


class CloningWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(CloningWidget, self).__init__(parent)

        self.detect_combo = QComboBox()
        self.detect_combo.addItems([self.tr('SIFT'), self.tr('BRISK')])
        self.detect_combo.setCurrentIndex(0)
        self.threshold_spin = QSpinBox()
        self.threshold_spin.setRange(1, 100)
        self.threshold_spin.setSuffix(self.tr('%'))
        self.distance_spin = QSpinBox()
        self.distance_spin.setRange(1, 100)
        self.distance_spin.setSuffix(self.tr('%'))
        self.distance_spin.setValue(10)
        self.cluster_spin = QSpinBox()
        self.cluster_spin.setRange(0, 20)
        self.cluster_spin.setValue(5)
        self.lines_check = QCheckBox(self.tr('Hide lines'))

        self.image = image
        self.viewer = ImageViewer(self.image, self.image)
        self.threshold = self.kpts = self.desc = self.matches = self.clusters = None
        self.gray = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        self.detect()

        self.detect_combo.currentIndexChanged.connect(self.detect)
        self.threshold_spin.valueChanged.connect(self.match)
        self.distance_spin.valueChanged.connect(self.cluster)
        self.cluster_spin.valueChanged.connect(self.cluster)
        self.lines_check.stateChanged.connect(self.draw)

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr('Detector:')))
        top_layout.addWidget(self.detect_combo)
        top_layout.addWidget(QLabel(self.tr('Threshold:')))
        top_layout.addWidget(self.threshold_spin)
        top_layout.addWidget(QLabel(self.tr('Distance:')))
        top_layout.addWidget(self.distance_spin)
        top_layout.addWidget(QLabel(self.tr('Cluster:')))
        top_layout.addWidget(self.cluster_spin)
        top_layout.addWidget(self.lines_check)
        top_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def detect(self):
        index = self.detect_combo.currentIndex()
        if index == 0:
            detector = cv.xfeatures2d.SIFT_create()
            self.threshold_spin.blockSignals(True)
            self.threshold_spin.setValue(40)
            self.threshold_spin.blockSignals(False)
        elif index == 1:
            detector = cv.BRISK_create()
            self.threshold_spin.blockSignals(True)
            self.threshold_spin.setValue(10)
            self.threshold_spin.blockSignals(False)
        else:
            return
        self.kpts, self.desc = detector.detectAndCompute(self.gray, None)
        self.match()

    def match(self):
        index = self.detect_combo.currentIndex()
        norm = cv.NORM_L2 if index == 0 else cv.NORM_HAMMING
        matcher = cv.BFMatcher_create(norm, True)
        self.threshold = self.threshold_spin.value() / 100 * 255
        self.matches = matcher.radiusMatch(self.desc, self.desc, self.threshold)
        self.matches = [item for sublist in self.matches for item in sublist]
        matches2 = []
        for m in self.matches:
            if m.queryIdx != m.trainIdx:
                matches2.append(m)
        self.matches = matches2
        self.cluster()

    def cluster(self):
        self.clusters = []
        total = len(self.matches)
        min_dist = (self.distance_spin.value() / 100) * (np.min(self.gray.shape) / 2)
        min_size = self.cluster_spin.value()
        for i in range(total):
            match0 = self.matches[i]
            query0 = match0.queryIdx
            train0 = match0.trainIdx
            group = [match0]
            a0 = np.array(self.kpts[query0].pt)
            b0 = np.array(self.kpts[train0].pt)
            if np.linalg.norm(a0 - b0) < min_dist:
                continue
            for j in range(i + 1, total):
                match1 = self.matches[j]
                query1 = match1.queryIdx
                train1 = match1.trainIdx
                if query1 == train0 and train1 == query0:
                    continue
                a1 = np.array(self.kpts[query1].pt)
                b1 = np.array(self.kpts[train1].pt)
                if np.linalg.norm(a1 - b1) < min_dist:
                    continue
                daa = np.linalg.norm(a0 - a1)
                dbb = np.linalg.norm(b0 - b1)
                dab = np.linalg.norm(a0 - b1)
                dba = np.linalg.norm(b0 - a1)
                smallest = np.partition(np.array([daa, dbb, dab, dba]), 1)[:2]
                if np.all(np.logical_and(smallest > 0, smallest < min_dist)):
                    for g in group:
                        if g.queryIdx == train1 and g.trainIdx == query1:
                            break
                    else:
                        group.append(match1)
            if len(group) >= min_size:
                self.clusters.append(group)
        self.draw()

    def draw(self):
        output = np.copy(self.image)
        hsv = np.array([[[0, 0, 255]]])
        lines = self.lines_check.isChecked()
        for i, c in enumerate(self.clusters):
            for m in c:
                ka = self.kpts[m.queryIdx]
                sa = int(np.round(ka.size))
                a = tuple(map(int, ka.pt))
                kb = self.kpts[m.trainIdx]
                sb = int(np.round(kb.size))
                b = tuple(map(int, kb.pt))
                if np.linalg.norm(a) > np.linalg.norm(b):
                    a, b = b, a
                hsv[0, 0, 0] = 90 / np.pi * (np.pi + np.arctan2(b[1] - a[1], b[0] - a[0]))
                hsv[0, 0, 1] = m.distance / self.threshold * 128 + 128
                rgb = cv.cvtColor(hsv.astype(np.uint8), cv.COLOR_HSV2BGR)
                rgb = tuple([int(x) for x in rgb[0, 0]])
                cv.circle(output, a, sa, rgb, 1, cv.LINE_AA)
                cv.circle(output, b, sb, rgb, 1, cv.LINE_AA)
                if not lines:
                    cv.line(output, a, b, rgb, 1, cv.LINE_AA)
        self.viewer.update_processed(output)
