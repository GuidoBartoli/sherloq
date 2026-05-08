from itertools import compress
from os.path import splitext
from time import time

from sklearn.cluster import DBSCAN
import cv2 as cv
import numpy as np
from PySide6.QtCore import QCoreApplication, Qt
from PySide6.QtWidgets import (
    QToolButton,
    QMessageBox,
    QCheckBox,
    QComboBox,
    QLabel,
    QHBoxLayout,
    QVBoxLayout,
    QSlider,
    QFrame,
)

from kneed import KneeLocator
from sklearn.neighbors import NearestNeighbors

from tools import ToolWidget
from utility import elapsed_time, modify_font, load_image
from viewer import ImageViewer


class CopyMoveWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(CopyMoveWidget, self).__init__(parent)

        # detector select
        self.detector_combo = QComboBox()
        self.detector_combo.addItems(
            [self.tr("BRISK"), self.tr("ORB"), self.tr("AKAZE"), self.tr("SIFT")]
        )
        self.detector_combo.setCurrentIndex(3)
        self.detector_combo.setToolTip(
            self.tr("Algorithm used for localization and description")
        )

        # response threshold
        self.response_slider = QSlider(Qt.Orientation.Horizontal)
        self.response_slider.setMinimum(0)
        self.response_slider.setMaximum(100)
        self.response_slider.setValue(100)
        self.response_slider.setToolTip(
            self.tr("Maximum keypoint response to perform matching")
        )

        # match threshold
        self.match_slider = QSlider(Qt.Orientation.Horizontal)
        self.match_slider.setMinimum(0)
        self.match_slider.setMaximum(100)
        self.match_slider.setValue(60)
        self.match_slider.setToolTip(self.tr("Ratio for controls match uniqueness"))

        # distance threshold (ratio * min(image.size))
        self.distance_slider = QSlider(Qt.Orientation.Horizontal)
        self.distance_slider.setMinimum(0)
        self.distance_slider.setMaximum(100)
        self.distance_slider.setValue(5)
        self.distance_slider.setToolTip(
            self.tr("Minimum physical distance between matched points")
        )

        # show keypoints
        self.kpts_check = QCheckBox(self.tr("Show keypoints"))
        self.kpts_check.setToolTip(self.tr("Show keypoint coverage"))

        # hide lines matches
        self.nolines_check = QCheckBox(self.tr("Hide lines"))
        self.nolines_check.setToolTip(self.tr("Disable match line drawing"))

        # hide region of tampered prediction
        self.noregion_check = QCheckBox(self.tr("Hide region area"))
        self.noregion_check.setToolTip(self.tr("Hide region area drawing"))

        # process button
        self.process_button = QToolButton()
        self.process_button.setText(self.tr("Process"))
        self.process_button.setToolTip(self.tr("Perform automatic detection"))
        modify_font(self.process_button, bold=True)

        self.status_label = QLabel()

        # load mask image
        self.mask_label = QLabel()
        self.mask_button = QToolButton()
        self.mask_button.setText(self.tr("Load mask..."))
        self.mask_button.setToolTip(self.tr("Load an image to be used as mask"))

        # onoff mask button
        self.onoff_button = QToolButton()
        self.onoff_button.setText(self.tr("OFF"))
        self.onoff_button.setCheckable(True)
        self.onoff_button.setToolTip(self.tr("Toggle keypoint detection mask"))

        self.image = image
        self.image = self.resize_image(self.image)

        self.viewer = ImageViewer(self.image, self.image)

        self.gray = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        self.total = self.kpts = self.desc = self.matches = self.mask = None
        self.canceled = False

        self.detector_combo.currentIndexChanged.connect(self.update_detector)
        self.response_slider.valueChanged.connect(self.update_detector)
        self.match_slider.valueChanged.connect(self.update_matching)
        self.distance_slider.valueChanged.connect(self.update_matching)
        self.nolines_check.stateChanged.connect(self.process)
        self.kpts_check.stateChanged.connect(self.process)
        self.process_button.clicked.connect(self.process)
        self.noregion_check.stateChanged.connect(self.process)
        self.mask_button.clicked.connect(self.load_mask)
        self.onoff_button.toggled.connect(self.toggle_mask)
        self.onoff_button.setEnabled(False)

        # response slider
        self.response_value_label = QLabel(f"{self.response_slider.value()}%")
        self.response_slider.valueChanged.connect(
            lambda value: self.response_value_label.setText(f"{value}%")
        )
        response_layout = QHBoxLayout()
        response_layout.addWidget(QLabel(self.tr("Response:")))
        response_layout.addWidget(self.response_slider)
        response_layout.addWidget(self.response_value_label)

        # match slider
        self.match_value_label = QLabel(f"{self.match_slider.value()}%")
        self.match_slider.valueChanged.connect(
            lambda value: self.match_value_label.setText(f"{value}%")
        )
        match_layout = QHBoxLayout()
        match_layout.addWidget(QLabel(self.tr("Match:")))
        match_layout.addWidget(self.match_slider)
        match_layout.addWidget(self.match_value_label)

        # distance slider
        self.distance_value_label = QLabel(f"{self.distance_slider.value()}%")
        self.distance_slider.valueChanged.connect(
            lambda value: self.distance_value_label.setText(f"{value}%")
        )
        distance_layout = QHBoxLayout()
        distance_layout.addWidget(QLabel(self.tr("Distance:")))
        distance_layout.addWidget(self.distance_slider)
        distance_layout.addWidget(self.distance_value_label)

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr("Detector:")))
        top_layout.addWidget(self.detector_combo)
        top_layout.addWidget(self.create_separator())
        top_layout.addLayout(response_layout)
        top_layout.addWidget(self.create_separator())
        top_layout.addLayout(match_layout)
        top_layout.addWidget(self.create_separator())
        top_layout.addLayout(distance_layout)

        middle_layout = QHBoxLayout()
        middle_layout.addStretch()
        middle_layout.addWidget(self.nolines_check)
        middle_layout.addWidget(self.create_separator())
        middle_layout.addWidget(self.noregion_check)
        middle_layout.addWidget(self.create_separator())
        middle_layout.addWidget(self.kpts_check)
        middle_layout.addStretch()

        bottom_layout = QHBoxLayout()
        bottom_layout.addWidget(self.process_button)
        bottom_layout.addWidget(self.status_label)
        bottom_layout.addStretch()
        bottom_layout.addWidget(self.mask_button)
        bottom_layout.addWidget(self.onoff_button)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addLayout(middle_layout)
        main_layout.addLayout(bottom_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def create_separator(self):
        # create separator
        separator = QFrame()
        separator.setFrameShape(QFrame.VLine)
        separator.setFrameShadow(QFrame.Sunken)
        return separator

    def toggle_mask(self, checked):
        self.onoff_button.setText("ON" if checked else "OFF")
        if checked:
            self.viewer.update_processed(
                cv.merge([c * self.mask for c in cv.split(self.image)])
            )
        else:
            self.viewer.update_processed(self.image)
        self.update_detector()

    def load_mask(self):
        filename, basename, mask = load_image(self)
        mask = self.resize_image(mask)  # resize mask
        if filename is None:
            return
        if self.image.shape[:-1] != mask.shape[:-1]:
            QMessageBox.critical(
                self,
                self.tr("Error"),
                self.tr("Both image and mask must have the same size!"),
            )
            return
        _, self.mask = cv.threshold(
            cv.cvtColor(mask, cv.COLOR_BGR2GRAY), 0, 1, cv.THRESH_BINARY
        )
        self.onoff_button.setEnabled(True)
        self.onoff_button.setChecked(True)
        self.mask_button.setText(f'"{splitext(basename)[0]}"')
        self.mask_button.setToolTip(self.tr("Current detection mask image"))

    def update_detector(self):
        # reset all
        self.total = self.kpts = self.desc = self.matches = None
        self.status_label.setText("")
        self.process_button.setEnabled(True)

    def update_matching(self):
        self.matches = None
        self.process_button.setEnabled(True)

    def cancel(self):
        self.canceled = True
        self.status_label.setText(self.tr("Processing interrupted!"))
        modify_font(self.status_label, bold=False, italic=False)

    def resize_image(self, img, max_dim=1920):
        # for large images, resize to speed up and reduce memory usage
        h, w = img.shape[:2]
        if max(h, w) > max_dim:
            self.status_label.setText(self.tr("Resizing image for processing"))
            scale = max_dim / max(h, w)
            return cv.resize(img, (int(w * scale), int(h * scale)))
        return img

    def _get_valid_pairs(self, dp, k_range=range(2, 11)):
        # get valid pairs of (eps, minPts) for final DBSCAN
        valid_pairs = []
        for k in k_range:
            if k > len(dp):
                break

            nn = NearestNeighbors(n_neighbors=k)
            nn.fit(dp)
            distances, _ = nn.kneighbors(dp)
            k_distances = np.sort(distances[:, -1])

            kneedle = KneeLocator(
                range(len(k_distances)),
                k_distances,
                S=1.0,
                curve="convex",
                direction="increasing",
            )

            if kneedle.knee is None:
                continue

            eps_candidate = float(k_distances[kneedle.knee])

            if eps_candidate <= 0.0:
                continue

            db = DBSCAN(eps=eps_candidate, min_samples=k)
            labels = db.fit_predict(dp)

            unique_labels = set(labels)
            n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)

            if 2 <= n_clusters <= 10:
                valid_pairs.append((k, eps_candidate))

        return valid_pairs

    def process(self):
        start = time()
        self.canceled = False
        self.status_label.setText(self.tr("Processing, please wait..."))

        algorithm = self.detector_combo.currentIndex()
        response = 100 - self.response_slider.value()
        match_ratio = self.match_slider.value() / 100.0
        distance_threshold = self.distance_slider.value() / 100.0

        modify_font(self.status_label, bold=False, italic=True)
        QCoreApplication.processEvents()

        # detector
        norm_type = cv.NORM_HAMMING
        desc_dtype = np.uint8
        if self.kpts is None:
            if algorithm == 0:
                detector = cv.BRISK_create(thresh=80)
            elif algorithm == 1:
                detector = cv.ORB_create(nfeatures=17000)
            elif algorithm == 2:
                detector = cv.AKAZE_create(threshold=0.0035)
            elif algorithm == 3:
                detector = cv.SIFT_create(nfeatures=17000)
                norm_type = cv.NORM_L2
                desc_dtype = np.float32
            else:
                return

            mask = self.mask if self.onoff_button.isChecked() else None

            # extract keypoints from gray image
            kp1, desc1 = detector.detectAndCompute(self.gray, mask)

            # extract keypoints from scaled scaled image (scale factor = 2.0)
            h, w = self.gray.shape[:2]
            scaled_gray = cv.resize(self.gray, (int(w * 2.0), int(h * 2.0)))
            scaled_kp2, desc2 = detector.detectAndCompute(scaled_gray, None)

            kp2 = []
            if scaled_kp2 is not None:
                for kp in scaled_kp2:
                    tmp_kp = cv.KeyPoint(
                        kp.pt[0] / 2.0,
                        kp.pt[1] / 2.0,
                        kp.size / 2.0,
                        kp.angle,
                        kp.response,
                        kp.octave,
                        kp.class_id,
                    )
                    kp2.append(tmp_kp)

            # combine keypoints and descriptors
            self.kpts = list(kp1) + kp2
            self.desc = np.vstack([desc1, desc2])

            if self.kpts is None or len(self.kpts) == 0:
                self.status_label.setText(self.tr("No keypoints found."))
                self.process_button.setEnabled(True)
                return

            self.total = len(self.kpts)

            responses = np.array([k.response for k in self.kpts])
            strongest = (
                cv.normalize(responses, None, 0, 100, cv.NORM_MINMAX) >= response
            ).flatten()
            self.kpts = list(compress(self.kpts, strongest))
            self.desc = self.desc[strongest]

        self.desc = self.desc.astype(desc_dtype)

        # second point matching
        if self.matches is None:
            matcher = cv.BFMatcher(normType=norm_type)
            self.matches = matcher.knnMatch(self.desc, self.desc, k=3)

            if self.matches is None:
                self.status_label.setText(
                    self.tr("No keypoint match found with current settings")
                )
                modify_font(self.status_label, italic=False, bold=True)
                return

            good_matches = []
            for m in self.matches:
                if len(m) < 3:
                    continue

                best_match = m[1]
                second_match = m[2]

                if best_match.distance < match_ratio * second_match.distance:
                    p1 = np.array(self.kpts[best_match.queryIdx].pt)
                    p2 = np.array(self.kpts[best_match.trainIdx].pt)
                    if np.linalg.norm(p1 - p2) > distance_threshold * min(
                        self.gray.shape[:2]
                    ):
                        good_matches.append((best_match.queryIdx, best_match.trainIdx))

            self.matches = good_matches

        # double filtering
        # 1. RANSAC homography
        if len(self.matches) >= 4:
            src_pts = np.float32([self.kpts[i].pt for i, _ in self.matches]).reshape(
                -1, 1, 2
            )
            dst_pts = np.float32([self.kpts[j].pt for _, j in self.matches]).reshape(
                -1, 1, 2
            )

            H, maskR = cv.findHomography(src_pts, dst_pts, cv.RANSAC, 50.0)
            maskR = maskR.ravel().astype(bool)

            inlier_matches = list(compress(self.matches, maskR))
        else:
            self.status_label.setText(self.tr("No valid matches found."))
            self.process_button.setEnabled(True)
            return

        dp_indices = np.unique(inlier_matches)
        dp_pts = np.array([self.kpts[i].pt for i in dp_indices])

        # 2. KANN-DBSCAN clustering
        valid_pairs = self._get_valid_pairs(dp_pts)

        if valid_pairs:
            minPts = [m for m, _ in valid_pairs]
            eps = [e for _, e in valid_pairs]

            # take the mean for final (eps, minPts)
            minPts_ada = int(np.mean(minPts))
            eps_ada = float(np.mean(eps))

            if eps_ada <= 0.0:
                eps_ada = 1e-5

            db = DBSCAN(eps=eps_ada, min_samples=minPts_ada).fit(dp_pts)
            labels = db.labels_
        else:
            labels = [-1] * len(dp_pts)

        # visualize results
        output = np.copy(self.image)
        overlay = output.copy()

        hide_regions = self.noregion_check.isChecked()
        valid_area = []
        if not hide_regions:
            for label in set(labels):
                if label == -1:
                    continue

                cluster_pts = np.array(
                    [dp_pts[i] for i in range(len(dp_pts)) if labels[i] == label]
                )

                # minium number of points to form a valid region is 3
                if len(cluster_pts) < 3:
                    continue

                hull = cv.convexHull(np.int32(cluster_pts))

                # convex hull is collinear
                if cv.contourArea(hull) == 0:
                    continue

                valid_area.append(hull)
                cv.fillConvexPoly(overlay, np.int32(hull), (0, 255, 0))
                cv.polylines(output, [np.int32(hull)], True, (0, 200, 0), 2, cv.LINE_4)

        alpha = 0.4
        cv.addWeighted(overlay, alpha, output, 1 - alpha, 0, output)

        nolines = self.nolines_check.isChecked()
        show_kpts = self.kpts_check.isChecked()

        if show_kpts:
            for kpt in self.kpts:
                cv.circle(output, (int(kpt.pt[0]), int(kpt.pt[1])), 2, (250, 227, 72))

        if not nolines:
            for i in range(len(inlier_matches)):
                idx1, idx2 = inlier_matches[i]
                p1 = (int(self.kpts[idx1].pt[0]), int(self.kpts[idx1].pt[1]))
                p2 = (int(self.kpts[idx2].pt[0]), int(self.kpts[idx2].pt[1]))
                cv.line(output, p1, p2, (255, 255, 0), 1, cv.LINE_4)

        self.viewer.update_processed(output)
        self.process_button.setEnabled(False)
        modify_font(self.status_label, italic=False, bold=True)
        self.status_label.setText(
            self.tr(
                f"Keypoints: {self.total} --> Inliers: {len(inlier_matches)} --> Regions: {len(valid_area)}"
            )
        )
        self.info_message.emit(self.tr(f"Copy-Move Forgery = {elapsed_time(start)}"))
