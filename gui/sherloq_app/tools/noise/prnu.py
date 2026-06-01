import os

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QPushButton,
    QSizePolicy,
    QStackedWidget,
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
        self.hdf5_path = os.path.dirname(
            f"{os.path.abspath(__file__)}/fingerprints/prnu.h5"
        )

        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(8, 8, 8, 8)
        main_layout.setSpacing(6)

        hdf5_layout = QHBoxLayout()
        self.hdf5_label = QLabel(self._short(self.hdf5_path))
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

        self.stack.addWidget(found_widget)
        self.stack.addWidget(not_found_widget)

        main_layout.addLayout(hdf5_layout)
        main_layout.addWidget(self._hline())
        main_layout.addWidget(self.stack)
        main_layout.addWidget(self.progress)
        main_layout.addWidget(self.status_label)

    def _browse_hdf5(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select PRNU Fingerprints File",
            os.path.dirname(self.hdf5_path),
            "HDF5 files (*.h5)",
        )

        if path:
            self.hdf5_path = path
            self.hdf5_label.setText(self._short(path))
            self.hdf5_label.setToolTip(path)
            self._check_hdf5()

    def _check_hdf5(self):
        if os.path.exists(self.hdf5_path):
            try:
                with h5py.File(self.hdf5_path, "r") as f:
                    n = len(f.keys())
                self.found_label.setText(
                    f"{n} camera fingerprint(s) ready\n{self._short(self.hdf5_path)}"
                )
            except:
                self.found_label.setText("Fingerprints file found")
            self.stack.setCurrentIndex(0)
            self.status_label.setText("Fingerprints loaded. Press Identify to run.")
        else:
            self.stack.setCurrentIndex(1)

    def _identify(self):
        return

    @staticmethod
    def _short(path: str, n: int = 55) -> str:
        """
        Shorten a file path for display purposes.
        """
        return ("…" + path[-(n - 1) :]) if len(path) > n else path

    @staticmethod
    def _hline():
        """
        Create and return a horizontal separator line widget.
        """
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)
        line.setFrameShadow(QFrame.Shadow.Sunken)
        return line
