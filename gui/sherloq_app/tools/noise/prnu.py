import os
import hashlib
from pathlib import Path

from PySide6.QtCore import QThread, Qt, Signal
from PySide6.QtGui import QColor, QFont
from PySide6.QtWidgets import (
    QAbstractItemView,
    QFileDialog,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)
from collections import defaultdict
from scipy.signal import wiener
import numpy as np
import h5py
import cv2

from gui.sherloq_app.ui.tools import ToolWidget

WIENER_SIZE = 3

def parse_camera_label(filename: str)-> str:
    stem = Path(filename).stem
    parts = stem.split("_")
    return "_".join(parts[:-1]) if len(parts) >= 3 else stem

def scan_dataset(path:str)-> dict:
    """Only cameras with at least 2 images are included in the result,
    since a single image is not enough to build a reliable PRNU fingerprint.
    """

    camera_imgs = defaultdict(list)
    ext = {".jpg", ".jpeg", ".JPG", ".JPEG"}
    for dirpath, _, files in os.walk(path):
        for f in files:
            if Path(f).suffix in ext:
                label = parse_camera_label(f)
                camera_imgs[label].append(os.path.join(dirpath, f))
    return {k: v for k, v in camera_imgs.items() if len(v) >= 2}

def extract_residual(image: np.ndarray) -> np.ndarray:
    return image - wiener(image, mysize=WIENER_SIZE)

def load_image_gray(path: str)-> np.ndarray:
    img = cv2.imread(path, cv2.IMREAD_COLOR)
    if img is None:
        raise IOError(f"Cannot read: {path}")
    return cv2.cvtColor(img, cv2.COLOR_BGR2GRAY).astype(np.float64) / 255.0

def get_folder_hash(paths: list) -> str:
    entries = [f"{os.path.basename(p)}:{os.path.getsize(p)}" for p in sorted(paths)]
    return hashlib.md5("\n".join(entries).encode()).hexdigest()

def ncc(a: np.ndarray, b: np.ndarray) -> float:
    h = min(a.shape[0], b.shape[0])
    w = min(a.shape[1], b.shape[1])
    a, b = a[:h, :w], b[:h, :w]
    a = a - a.mean()
    b = b - b.mean()
    denom = np.sqrt((a**2).sum() * (b**2).sum())
    return float(np.sum(a * b) / denom) if denom > 1e-10 else 0.0

class PrnuWorker(QThread):
    """Runs PRNU work off the UI thread."""
    progress = Signal(int, str)
    finished = Signal(list)
    error = Signal(str)
    cancelled = Signal() 

    def __init__(
        self, image: np.ndarray, hdf5_path: str, dataset_dir: str | None = None
    ) -> None:
        super().__init__()
        self.image = image
        self.hdf5_path = hdf5_path
        self.dataset_dir = dataset_dir
        self.cancel_flag = False

    def request_cancel(self):
        """Safe to call from any thread. Stops at the next image boundary."""
        self.cancel_flag = True

    def run(self):
        try:
            fingerprints = {}
            if os.path.exists(self.hdf5_path):
                self.progress.emit(5, "Loading fingerprints from HDF5")
                with h5py.File(self.hdf5_path, "r") as f:
                    cams = list(f.keys())
                    if not cams:
                        self.error.emit("HDF5 file is empty")
                        return
                    n = len(cams)
                    for i, cam in enumerate(cams):
                        self.progress.emit(5 + int(30 * (i + 1) / n), f"Loading: {cam}")
                        fingerprints[cam] = f[cam]["fingerprint"][:] 

            else:
                if not self.dataset_dir:
                    self.error.emit("No dataset folder provided.")
                    return

                self.progress.emit(2, "Scanning dataset folder")
                camera_images = scan_dataset(self.dataset_dir)
                if not camera_images:
                    self.error.emit(
                        f"No valid images found in:\n{self.dataset_dir}\n\n"
                    )
                    return

                n_cams = len(camera_images)
                self.progress.emit(
                    5, f"Found {n_cams} camera(s). Building fingerprints ..."
                )

                for cam_idx, (cam, paths) in enumerate(camera_images.items()):
                    cam_start = 5 + int(60 * cam_idx / n_cams)
                    cam_end = 5 + int(60 * (cam_idx + 1) / n_cams)
                    mean_fp = None
                    n_used = 0

                    for img_idx, path in enumerate(paths):
                        if self.cancel_flag:
                            self.cancelled.emit()
                            return

                        pct = cam_start + int(
                            (cam_end - cam_start) * (img_idx + 1) / len(paths)
                        )
                        self.progress.emit(
                            pct,
                            f"[{cam_idx+1}/{n_cams}] {cam}  "
                            f"image {img_idx+1}/{len(paths)} ...",
                        )
                        try:
                            residual = extract_residual(load_image_gray(path))
                            if mean_fp is None:
                                mean_fp = residual.copy()
                                n_used = 1
                            else:
                                h = min(residual.shape[0], mean_fp.shape[0])
                                w = min(residual.shape[1], mean_fp.shape[1])
                                mean_fp = mean_fp[:h, :w]
                                n_used += 1
                                mean_fp += (residual[:h, :w] - mean_fp) / n_used
                        except Exception:
                            continue

                    if self.cancel_flag:
                        self.cancelled.emit()
                        return

                    if mean_fp is not None:
                        fingerprints[cam] = mean_fp
                        self.progress.emit(
                            cam_end, f"[{cam_idx+1}/{n_cams}] {cam} - saving to HDF5"
                        )
                        with h5py.File(self.hdf5_path, "a") as f:
                            if cam in f:
                                del f[cam]
                            grp = f.create_group(cam)
                            grp.create_dataset(
                                "fingerprint",
                                data=mean_fp,
                                dtype=np.float64,
                                compression="gzip",
                                compression_opts=4,
                            )
                            grp.attrs["n_images"] = len(paths)
                            grp.attrs["n_used"] = n_used
                            grp.attrs["shape"] = list(mean_fp.shape)
                            grp.attrs["folder_hash"] = get_folder_hash(paths)
                        self.progress.emit(
                            cam_end,
                            f"[{cam_idx+1}/{n_cams}] {cam} - saved"
                            f"({n_used}/{len(paths)} images)",
                        )

                self.progress.emit(
                    65, f"All {len(fingerprints)} fingerprint(s) saved to HDF5"
                )

            if not fingerprints:
                self.error.emit("No fingerprints available.")
                return

            self.progress.emit(68, "Extracting PRNU residual from query image")
            gray = (
                cv2.cvtColor(self.image, cv2.COLOR_BGR2GRAY).astype(np.float64) / 255.0
            )
            residual = extract_residual(gray)

            results = []
            n = len(fingerprints)
            for i, (cam, fp) in enumerate(fingerprints.items()):
                self.progress.emit(70 + int(28 * (i + 1) / n), f"Matching: {cam}")
                results.append((cam, ncc(residual, fp)))

            results.sort(key=lambda x: x[1], reverse=True)
            self.progress.emit(100, "Identification complete.")
            self.finished.emit(results)

        except Exception as e:
            self.error.emit(str(e))


class PrnuWidget(ToolWidget):
    def __init__(self, filename: str, image: np.ndarray, parent=None):
        super().__init__(parent)
        self.filename = filename
        self.image = image
        self.hdf5_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "fingerprints", "prnu.h5"
        )
        self.dataset_dir = None
        self.worker = None
        self._on_start()
        self._check_hdf5()

    def _on_start(self):
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(8, 8, 8, 8)
        main_layout.setSpacing(6)

        hdf5_layout = QHBoxLayout()
        self.hdf5_label = QLabel(self.hdf5_path)
        self.hdf5_label.setToolTip(self.hdf5_path)
        self.hdf5_label.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        browse_h5_button = QPushButton("Browse")
        browse_h5_button.setFixedWidth(90)
        browse_h5_button.clicked.connect(self._browse_hdf5)
        hdf5_layout.addWidget(QLabel("Fingerprints (.h5):"))
        hdf5_layout.addWidget(self.hdf5_label)
        hdf5_layout.addWidget(browse_h5_button)

        self.stack = QStackedWidget()

        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.status_label = QLabel("Ready.")

        # index 0 (h5 found)
        found_widget = QWidget()
        found_layout = QVBoxLayout(found_widget)
        found_layout.setContentsMargins(0, 0, 0, 0)
        self.found_label = QLabel()
        self.found_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.found_label.setStyleSheet("color: #1a7f1a; font-weight: bold;")
        self.identify_button = QPushButton("Identify Camera Source")
        self.identify_button.setFixedHeight(34)
        self.identify_button.clicked.connect(self._identify)

        found_layout.addWidget(self.found_label)
        found_layout.addWidget(self.identify_button)

        # index 1 (there is no h5)
        not_found_widget = QWidget()
        nf_layout = QVBoxLayout(not_found_widget)
        nf_layout.setContentsMargins(0, 0, 0, 0)
        default_label = QLabel(
            "No fingerprints database found.\n"
            "Select your camera image dataset folder to generate fingerprints, then identification runs automatically."
        )
        default_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        default_label.setStyleSheet("color: #b85c00;")

        dataset_layout = QHBoxLayout()
        self.dataset_label = QLabel("No folder selected")
        self.dataset_label.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        browse_dataset_button = QPushButton("Select Dataset")
        browse_dataset_button.setFixedWidth(130)
        browse_dataset_button.clicked.connect(self._browse_dataset)
        dataset_layout.addWidget(QLabel("Dataset folder:"))
        dataset_layout.addWidget(self.dataset_label)
        dataset_layout.addWidget(browse_dataset_button)

        gen_fp_layout = QHBoxLayout()
        self.gen_fp_button = QPushButton("Generate Fingerprints & Identify")
        self.gen_fp_button.setFixedHeight(34)
        self.gen_fp_button.setEnabled(False)
        self.gen_fp_button.clicked.connect(self._identify)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setFixedHeight(34)
        self.cancel_button.setFixedWidth(100)
        self.cancel_button.setEnabled(False)
        self.cancel_button.setStyleSheet("""
            QPushButton          { color: #aaaaaa; font-weight: bold;
                                   background-color: #e8e8e8;
                                   border: 1px solid #cccccc; border-radius: 4px; }
            QPushButton:enabled  { color: #c0392b; background-color: #fdecea;
                                   border: 1px solid #e74c3c; }
            QPushButton:enabled:hover   { background-color: #f5b7b1; }
            QPushButton:enabled:pressed { background-color: #e74c3c; color: white; }
        """)
        self.cancel_button.clicked.connect(self._cancel)
        gen_fp_layout.addWidget(self.gen_fp_button)
        gen_fp_layout.addWidget(self.cancel_button)

        nf_layout.addWidget(default_label)
        nf_layout.addLayout(dataset_layout)
        nf_layout.addLayout(gen_fp_layout)

        self.stack.addWidget(found_widget)
        self.stack.addWidget(not_found_widget)

        # result of identification
        res_box = QGroupBox("Identification Results")
        res_layout = QVBoxLayout(res_box)
        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(
            ["Rank", "Camera", "NCC Score", "Confidence"]
        )
        self.table.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeMode.Stretch
        )
        self.table.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.ResizeMode.ResizeToContents
        )
        self.table.horizontalHeader().setSectionResizeMode(
            3, QHeaderView.ResizeMode.ResizeToContents
        )
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setAlternatingRowColors(True)
        self.table.verticalHeader().setVisible(False)
        res_layout.addWidget(self.table)

        self.verdict = QLabel("")
        self.verdict.setAlignment(Qt.AlignmentFlag.AlignCenter)
        vf = QFont()
        vf.setPointSize(10)
        vf.setBold(True)
        self.verdict.setFont(vf)

        main_layout.addLayout(hdf5_layout)
        main_layout.addWidget(self.hline())
        main_layout.addWidget(self.stack)
        main_layout.addWidget(self.progress)
        main_layout.addWidget(self.status_label)
        main_layout.addWidget(res_box)
        main_layout.addWidget(self.verdict)

    def _browse_hdf5(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select PRNU Fingerprints File",
            os.path.dirname(self.hdf5_path),
            "HDF5 files (*.h5)",
        )

        if path:
            self.hdf5_path = path
            self.hdf5_label.setText(path)
            self.hdf5_label.setToolTip(path)
            self._check_hdf5()

    def _check_hdf5(self):
        if os.path.exists(self.hdf5_path):
            try:
                with h5py.File(self.hdf5_path, "r") as f:
                    n = len(f.keys())
                self.found_label.setText(
                    f"{n} camera fingerprint(s) ready\n{self.hdf5_path}"
                )
            except:
                self.found_label.setText("Fingerprints file found")
            self.stack.setCurrentIndex(0)
            self.status_label.setText("Fingerprints loaded. Press Identify to run.")
        else:
            self.stack.setCurrentIndex(1)

    def _identify(self):
        if self.image is None:
            QMessageBox.warning(self, "No Image", "No Image loaded in Sherloq")
            return
        if not os.path.exists(self.hdf5_path) and not self.dataset_dir:
            QMessageBox.warning(self, "Dataset required", "No fingerprints file and no dataset folder selected")
            return

        self._set_busy(True)
        self.table.setRowCount(0)
        self.verdict.setText("")
        self.progress.setValue(0)
        self.info_message.emit("PRNU identification is running ...")

        self.worker = PrnuWorker(self.image, self.hdf5_path, self.dataset_dir)
        self.worker.progress.connect(self._on_progress)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)
        self.worker.cancelled.connect(self._on_cancelled)
        self.worker.start()

    def _set_busy(self, busy: bool):
        """Toggle buttons correctly for running vs idle state"""
        self.identify_button.setEnabled(not busy)
        self.gen_fp_button.setEnabled(not busy and bool(self.dataset_dir))
        generating = busy and not os.path.exists(self.hdf5_path)
        self.cancel_button.setEnabled(generating)

    def _browse_dataset(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Select Camera Image Dataset Folder", os.path.expanduser("~")
        )
        if folder:
            self.dataset_dir = folder
            self.dataset_label.setText(folder)
            self.dataset_label.setToolTip(folder)
            self.gen_fp_button.setEnabled(True)

    def _cancel(self):
        """Request cancellation, worker stops at next image boundary."""
        if self.worker and self.worker.isRunning():
            self.cancel_button.setEnabled(False)
            self.status_label.setText("Cancelling... finishing current image...")
            self.worker.request_cancel()

    def _on_progress(self, pct, msg):
        self.progress.setValue(pct)
        self.status_label.setText(msg)

    def _on_finished(self, res):
        self._set_busy(False)
        self.progress.setValue(100)

        if self.stack.currentIndex() == 1  and os.path.exists(self.hdf5_path):
            self._check_hdf5()

        if not res:
            self.status_label.setText("No results.")
            return

        top = res[0][1]
        sec = res[1][1] if len(res) > 1 else 0.0
        margin = top - sec


        self.table.setRowCount(len(res))
        for row, (cam, score) in enumerate(res):
            conf = ("Strong" if margin > 0.01 else "Weak") if row == 0 else ""
            items = [
                QTableWidgetItem(str(row + 1)),
                QTableWidgetItem(cam),
                QTableWidgetItem(f"{score:.5f}"),
                QTableWidgetItem(conf),
            ]
            for item in items:
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
            if row == 0:
                for item in items:
                    item.setBackground(QColor(200, 255, 200))
            for col, item in enumerate(items):
                self.table.setItem(row, col, item)

        NCC_THRESHOLD = 0.005 
        if top < NCC_THRESHOLD:
            self.verdict.setText(f"Camera is not found in database")
            self.verdict.setStyleSheet("color: #8e44ad; font-weight: bold;")
            self.status_label.setText(
                f"The source camera may not be in the fingerprint database."
            )
            self.info_message.emit("PRNU: camera is not in database.")
            return

        self.verdict.setText(f"Predicted source camera:  {res[0][0]}   ")
        self.verdict.setStyleSheet(
            "color: #1a7f1a;" if margin > 0.01 else "color: #c47a00;"
        )
        self.status_label.setText(
            f"Matched {len(res)} camera(s)"
        )
        self.info_message.emit(f"PRNU: predicted camera ΓåÆ {res[0][0]}")

    def _on_cancelled(self):
        if self.worker:
            self.worker.wait()

        if os.path.exists(self.hdf5_path):
            try:
                os.remove(self.hdf5_path)
            except Exception:
                self.status_label.setText(f"Could not delete partial HDF5")
                self._set_busy(False)
                return

        self.progress.setValue(0)
        self.table.setRowCount(0)
        self.verdict.setText("")
        self._set_busy(False)
        self._check_hdf5()

    def _on_error(self, msg):
        self._set_busy(False)
        self.progress.setValue(0)
        self.status_label.setText(f"Error: {msg}")
        QMessageBox.critical(self, "PRNU Error", msg)
        self.info_message.emit(f"PRNU identification failed.")

    @staticmethod
    def hline():
        """
        Create and return a horizontal separator line widget.
        """
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        return line
