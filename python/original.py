from PySide2.QtWidgets import QWidget, QVBoxLayout

from viewer import ImageViewer


class OriginalWidget(QWidget):
    def __init__(self, filename, image, parent=None):
        super(OriginalWidget, self).__init__(parent)
        viewer = ImageViewer(image, None, None)
        layout = QVBoxLayout()
        layout.addWidget(viewer)
        self.setLayout(layout)
