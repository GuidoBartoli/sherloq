from time import time

import cv2 as cv
import numpy as np
from PySide2.QtCore import Qt
from PySide2.QtWidgets import (
    QPushButton,
    QSpinBox,
    QCheckBox,
    QComboBox,
    QLabel,
    QHBoxLayout,
    QVBoxLayout,
    QProgressDialog)

from tools import ToolWidget
from utility import elapsed_time
from viewer import ImageViewer


class CloningWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(CloningWidget, self).__init__(parent)

        self.detector_combo = QComboBox()
        self.detector_combo.addItems([self.tr('BRISK'), self.tr('ORB')])
        try:
            cv.xfeatures2d.SIFT_create()
            cv.xfeatures2d.SURF_create()
            self.detector_combo.addItem(self.tr('SIFT'))
            self.detector_combo.addItem(self.tr('SURF'))
        except cv.error:
            print('WARNING: Patented SIFT/SURF detectors not enabled in OpenCV package')
        self.detector_combo.setCurrentIndex(0)
        self.threshold_spin = QSpinBox()
        self.threshold_spin.setRange(1, 100)
        self.threshold_spin.setSuffix(self.tr('%'))
        self.threshold_spin.setValue(20)
        self.distance_spin = QSpinBox()
        self.distance_spin.setRange(1, 100)
        self.distance_spin.setSuffix(self.tr('%'))
        self.distance_spin.setValue(10)
        self.cluster_spin = QSpinBox()
        self.cluster_spin.setRange(0, 20)
        self.cluster_spin.setValue(4)
        self.nolines_check = QCheckBox(self.tr('Hide lines'))
        self.process_button = QPushButton(self.tr('Process'))

        self.image = image
        self.viewer = ImageViewer(self.image, self.image)
        self.gray = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        self.kpts = self.desc = self.matches = self.clusters = None
        self.canceled = False

        self.detector_combo.currentIndexChanged.connect(self.update_detector)
        self.threshold_spin.valueChanged.connect(self.update_threshold)
        self.distance_spin.valueChanged.connect(self.update_cluster)
        self.cluster_spin.valueChanged.connect(self.update_cluster)
        self.nolines_check.stateChanged.connect(self.process)
        self.process_button.clicked.connect(self.process)

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr('Detector:')))
        top_layout.addWidget(self.detector_combo)
        top_layout.addWidget(QLabel(self.tr('Threshold:')))
        top_layout.addWidget(self.threshold_spin)
        top_layout.addWidget(QLabel(self.tr('Distance:')))
        top_layout.addWidget(self.distance_spin)
        top_layout.addWidget(QLabel(self.tr('Cluster:')))
        top_layout.addWidget(self.cluster_spin)
        top_layout.addWidget(self.nolines_check)
        top_layout.addWidget(self.process_button)
        top_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def update_detector(self):
        self.kpts = self.desc = self.matches = self.clusters = None
        self.process_button.setEnabled(True)

    def update_threshold(self):
        self.matches = self.clusters = None
        self.process_button.setEnabled(True)

    def update_cluster(self):
        self.clusters = None
        self.process_button.setEnabled(True)

    def cancel(self):
        self.canceled = True

    def process(self):
        start = time()
        threshold = self.threshold_spin.value() * 2

        if self.kpts is None:
            index = self.detector_combo.currentIndex()
            if index == 0:
                detector = cv.BRISK_create()
            elif index == 1:
                detector = cv.ORB_create()
            elif index == 2:
                detector = cv.xfeatures2d.SURF_create(extended=True, upright=True)
            else:
                detector = cv.xfeatures2d.SIFT_create()
            self.kpts, self.desc = detector.detectAndCompute(self.gray, None)

        if self.matches is None:
            index = self.detector_combo.currentIndex()
            norm = cv.NORM_HAMMING if index <= 1 else cv.NORM_L2
            matcher = cv.BFMatcher_create(norm, True)
            self.matches = matcher.radiusMatch(self.desc, self.desc, threshold)
            self.matches = [item for sublist in self.matches for item in sublist]
            self.matches = [m for m in self.matches if m.queryIdx != m.trainIdx]

        if self.clusters is None:
            self.clusters = []
            total = len(self.matches)
            min_dist = (self.distance_spin.value() / 100) * (np.min(self.gray.shape) / 2)
            min_size = self.cluster_spin.value()
            progress = QProgressDialog(self.tr('Clustering matches...'), self.tr('Cancel'), 0, total, self)
            progress.canceled.connect(self.cancel)
            progress.setWindowModality(Qt.WindowModal)
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
                    aa = np.linalg.norm(a0 - a1)
                    bb = np.linalg.norm(b0 - b1)
                    ab = np.linalg.norm(a0 - b1)
                    ba = np.linalg.norm(b0 - a1)
                    smallest = np.partition(np.array([aa, bb, ab, ba]), 1)[:2]
                    if np.all(np.logical_and(smallest > 0, smallest < min_dist)):
                        for g in group:
                            if g.queryIdx == train1 and g.trainIdx == query1:
                                break
                        else:
                            group.append(match1)
                if len(group) >= min_size:
                    self.clusters.append(group)
                progress.setValue(i)
                if self.canceled:
                    return
            progress.setValue(total)

        output = np.copy(self.image)
        hsv = np.zeros((1, 1, 3))
        nolines = self.nolines_check.isChecked()
        for c in self.clusters:
            for m in c:
                ka = self.kpts[m.queryIdx]
                pa = tuple(map(int, ka.pt))
                sa = int(np.round(ka.size))
                kb = self.kpts[m.trainIdx]
                pb = tuple(map(int, kb.pt))
                sb = int(np.round(kb.size))
                angle = np.arctan2(pb[1] - pa[1], pb[0] - pa[0])
                if angle < 0:
                    angle += np.pi
                hsv[0, 0, 0] = angle / np.pi * 180
                hsv[0, 0, 1] = 255
                hsv[0, 0, 2] = m.distance / threshold * 255
                rgb = cv.cvtColor(hsv.astype(np.uint8), cv.COLOR_HSV2BGR)
                rgb = tuple([int(x) for x in rgb[0, 0]])
                cv.circle(output, pa, sa, rgb, 1, cv.LINE_AA)
                cv.circle(output, pb, sb, rgb, 1, cv.LINE_AA)
                if not nolines:
                    cv.line(output, pa, pb, rgb, 1, cv.LINE_AA)

        self.viewer.update_processed(output)
        self.process_button.setEnabled(False)
        self.info_message.emit(self.tr('Copy-Move Forgery = {}'.format(elapsed_time(start))))
