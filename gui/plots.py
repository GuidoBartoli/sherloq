from PySide2.QtCharts.QtCharts import QScatterSeries
from PySide2.QtWidgets import (
    QWidget)


class PlotsWidget(QWidget):
    def __init__(self, image, parent=None):
        super(PlotsWidget, self).__init__(parent)
        series = QScatterSeries()
        series.setName("test")


