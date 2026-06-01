from PySide6.QtWidgets import QVBoxLayout
import numpy as np

from gui.sherloq_app.ui.tools import ToolWidget


class PrnuWidget(ToolWidget):
    def __init__(self, filename: str, image: np.ndarray, parent=None):
        super().__init__(parent)
        self.filename = filename
        self.image = image

        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(8, 8, 8, 8)
        main_layout.setSpacing(6)
