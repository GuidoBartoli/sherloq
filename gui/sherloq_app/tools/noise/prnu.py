import os

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QAbstractItemView,
    QFileDialog,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QStackedWidget,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)
import numpy as np
import h5py

from gui.sherloq_app.ui.tools import ToolWidget


class PrnuWidget(ToolWidget):
    def __init__(self, filename: str, image: np.ndarray, parent=None):
        super().__init__(parent)
        self.filename = filename
        self.image = image
        self.hdf5_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "fingerprints", "prnu.h5"
        )
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
        main_layout.addWidget(self._hline())
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
        return

    def _browse_dataset(self):
        return

    def _cancel(self):
        return

    @staticmethod
    def _hline():
        """
        Create and return a horizontal separator line widget.
        """
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        return line
