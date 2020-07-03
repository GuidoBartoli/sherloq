from itertools import compress
from time import time

import cv2 as cv
import numpy as np
from PySide2.QtCore import Qt, QCoreApplication
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
from utility import elapsed_time, modify_font
from viewer import ImageViewer


class CloningWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(CloningWidget, self).__init__(parent)

        self.detector_combo = QComboBox()
        self.detector_combo.addItems([self.tr('BRISK'), self.tr('ORB'), self.tr('AKAZE')])
        self.detector_combo.setCurrentIndex(0)
        self.detector_combo.setToolTip(self.tr('Algorithm used for localization and description'))
        self.response_spin = QSpinBox()
        self.response_spin.setRange(0, 100)
        self.response_spin.setSuffix(self.tr('%'))
        self.response_spin.setValue(90)
        self.response_spin.setToolTip(self.tr('Maximum keypoint response to perform matching'))
        self.matching_spin = QSpinBox()
        self.matching_spin.setRange(1, 100)
        self.matching_spin.setSuffix(self.tr('%'))
        self.matching_spin.setValue(20)
        self.matching_spin.setToolTip(self.tr('Maximum metric difference to accept matching'))
        self.distance_spin = QSpinBox()
        self.distance_spin.setRange(1, 100)
        self.distance_spin.setSuffix(self.tr('%'))
        self.distance_spin.setValue(15)
        self.distance_spin.setToolTip(self.tr('Maximum distance between matches in the same cluster'))
        self.cluster_spin = QSpinBox()
        self.cluster_spin.setRange(1, 20)
        self.cluster_spin.setValue(5)
        self.cluster_spin.setToolTip(self.tr('Minimum number of keypoints to create a new cluster'))
        self.nolines_check = QCheckBox(self.tr('Hide lines'))
        self.nolines_check.setToolTip(self.tr('Disable match line drawing'))
        self.process_button = QPushButton(self.tr('Process'))
        self.process_button.setToolTip(self.tr('Perform automatic detection'))
        self.status_label = QLabel(self.tr('[Press "Process" button to search for cloned regions]'))

        self.image = image
        self.viewer = ImageViewer(self.image, self.image)
        self.gray = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        self.keypoints = self.kpts = self.desc = self.matches = self.clusters = None
        self.canceled = False

        self.detector_combo.currentIndexChanged.connect(self.update_detector)
        self.response_spin.valueChanged.connect(self.update_detector)
        self.matching_spin.valueChanged.connect(self.update_matching)
        self.distance_spin.valueChanged.connect(self.update_cluster)
        self.cluster_spin.valueChanged.connect(self.update_cluster)
        self.nolines_check.stateChanged.connect(self.process)
        self.process_button.clicked.connect(self.process)

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr('Detector:')))
        top_layout.addWidget(self.detector_combo)
        top_layout.addWidget(QLabel(self.tr('Response:')))
        top_layout.addWidget(self.response_spin)
        top_layout.addWidget(QLabel(self.tr('Matching:')))
        top_layout.addWidget(self.matching_spin)
        top_layout.addWidget(QLabel(self.tr('Distance:')))
        top_layout.addWidget(self.distance_spin)
        top_layout.addWidget(QLabel(self.tr('Cluster:')))
        top_layout.addWidget(self.cluster_spin)
        top_layout.addWidget(self.nolines_check)
        top_layout.addWidget(self.process_button)
        top_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.status_label)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def update_detector(self):
        self.keypoints = self.kpts = self.desc = self.matches = self.clusters = None
        self.process_button.setEnabled(True)

    def update_matching(self):
        self.matches = self.clusters = None
        self.process_button.setEnabled(True)

    def update_cluster(self):
        self.clusters = None
        self.process_button.setEnabled(True)

    def cancel(self):
        self.canceled = True
        self.keypoints = self.kpts = self.desc = self.matches = self.clusters = None
        self.status_label.setText(self.tr('Processing interrupted!'))
        modify_font(self.status_label, bold=False, italic=False)

    def process(self):
        start = time()
        self.status_label.setText(self.tr('Processing, please wait...'))
        algorithm = self.detector_combo.currentIndex()
        response = 100 - self.response_spin.value()
        matching = self.matching_spin.value() / 100 * 255
        distance = self.distance_spin.value() / 100
        cluster = self.cluster_spin.value()
        modify_font(self.status_label, bold=False, italic=True)
        QCoreApplication.processEvents()

        if self.kpts is None:
            if algorithm == 0:
                detector = cv.BRISK_create()
            elif algorithm == 1:
                detector = cv.ORB_create()
            elif algorithm == 2:
                detector = cv.AKAZE_create()
            else:
                return
            self.kpts, self.desc = detector.detectAndCompute(self.gray, None)
            self.keypoints = len(self.kpts)
            responses = np.array([k.response for k in self.kpts])
            strongest = (cv.normalize(responses, None, 0, 100, cv.NORM_MINMAX) >= response).flatten()
            self.kpts = list(compress(self.kpts, strongest))
            self.desc = self.desc[strongest]

        if self.matches is None:
            matcher = cv.BFMatcher_create(cv.NORM_HAMMING, True)
            self.matches = matcher.radiusMatch(self.desc, self.desc, matching)
            if self.matches is None:
                self.status_label.setText(self.tr('No keypoint match found with current settings'))
                modify_font(self.status_label, italic=False, bold=True)
                return
            self.matches = [item for sublist in self.matches for item in sublist]
            self.matches = [m for m in self.matches if m.queryIdx != m.trainIdx]

        if self.clusters is None:
            self.clusters = []
            total = len(self.matches)
            min_dist = distance * np.min(self.gray.shape) / 2
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
                d0 = np.linalg.norm(a0 - b0)
                if d0 < min_dist:
                    continue
                for j in range(i + 1, total):
                    match1 = self.matches[j]
                    query1 = match1.queryIdx
                    train1 = match1.trainIdx
                    if query1 == train0 and train1 == query0:
                        continue
                    a1 = np.array(self.kpts[query1].pt)
                    b1 = np.array(self.kpts[train1].pt)
                    d1 = np.linalg.norm(a1 - b1)
                    if d1 < min_dist or np.abs(d0 - d1) > min_dist:
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
                if len(group) >= cluster:
                    self.clusters.append(group)
                progress.setValue(i)
                if self.canceled:
                    self.canceled = False
                    return
            progress.setValue(total)

        output = np.copy(self.image)
        hsv = np.zeros((1, 1, 3))
        nolines = self.nolines_check.isChecked()
        angles = []
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
                angles.append(angle)
                hsv[0, 0, 0] = angle / np.pi * 180
                hsv[0, 0, 1] = 255
                hsv[0, 0, 2] = m.distance / matching * 255
                rgb = cv.cvtColor(hsv.astype(np.uint8), cv.COLOR_HSV2BGR)
                rgb = tuple([int(x) for x in rgb[0, 0]])
                cv.circle(output, pa, sa, rgb, 1, cv.LINE_AA)
                cv.circle(output, pb, sb, rgb, 1, cv.LINE_AA)
                if not nolines:
                    cv.line(output, pa, pb, rgb, 1, cv.LINE_AA)

        regions = 0
        if angles:
            angles = np.reshape(np.array(angles, dtype=np.float32), (len(angles), 1))
            if np.std(angles) < 0.1:
                regions = 1
            else:
                criteria = (cv.TERM_CRITERIA_EPS + cv.TERM_CRITERIA_MAX_ITER, 10, 1.0)
                attempts = 10
                flags = cv.KMEANS_PP_CENTERS
                compact = [cv.kmeans(angles, k, None, criteria, attempts, flags)[0] for k in range(1, 11)]
                compact = cv.normalize(np.array(compact), None, 0, 1, cv.NORM_MINMAX)
                regions = np.argmax(compact < 0.005) + 1
        self.viewer.update_processed(output)
        self.process_button.setEnabled(False)
        modify_font(self.status_label, italic=False, bold=True)
        self.status_label.setText(
            self.tr('Keypoints: {} --> Filtered: {} --> Matches: {} --> Clusters: {} --> Regions: {}'.format(
                self.keypoints, len(self.kpts), len(self.matches), len(self.clusters), regions)))
        self.info_message.emit(self.tr('Copy-Move Forgery = {}'.format(elapsed_time(start))))
