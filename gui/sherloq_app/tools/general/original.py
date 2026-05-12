from PySide6.QtWidgets import QVBoxLayout

from gui.sherloq_app.ui.tools import ToolWidget
from gui.sherloq_app.ui.viewer import ImageViewer


class OriginalWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(OriginalWidget, self).__init__(parent)
        viewer = ImageViewer(image, None, None)
        layout = QVBoxLayout()
        layout.addWidget(viewer)
        self.setLayout(layout)
