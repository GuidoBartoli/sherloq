import cv2 as cv
from PySide2.QtWidgets import (
    QWidget,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QSpinBox,
    QSlider,
    QLabel)
from viewer import ImageViewer


class PcaWidget(QWidget):
    def __init__(self, image, parent=None):
        super(PcaWidget, self).__init__(parent)

