import cv2 as cv
from PySide2.QtWidgets import (
    QWidget,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QHBoxLayout,
    QSpinBox,
    QSlider,
    QLabel)
from PySide2.QtCore import Qt
from PySide2.QtCharts.QtCharts import QScatterSeries

from utility import compress_jpeg
from viewer import ImageViewer


class PlotsWidget(QWidget):
    def __init__(self, image, parent=None):
        super(PlotsWidget, self).__init__(parent)
        series = QScatterSeries()
        series.setName("test")


