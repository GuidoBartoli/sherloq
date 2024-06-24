import cv2 as cv
from PySide6.QtCharts import QLineSeries, QChart, QChartView
from PySide6.QtCore import Qt
from PySide6.QtGui import QPainter
from PySide6.QtWidgets import QVBoxLayout, QProgressDialog

from jpeg import compress_jpg
from tools import ToolWidget


class MultipleWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(ToolWidget, self).__init__(parent)

        max_q = 101
        progress = QProgressDialog(
            self.tr("Computing residuals..."), None, 0, max_q, self
        )
        progress.setWindowModality(Qt.WindowModal)
        loss_series = QLineSeries()
        gray = cv.cvtColor(image, cv.COLOR_BGR2GRAY)
        for q in range(max_q):
            loss = cv.mean(cv.absdiff(compress_jpg(gray, q, color=False), gray))
            loss_series.append(q, loss[0])
            progress.setValue(q)
        progress.setValue(max_q)

        loss_chart = QChart()
        loss_chart.legend().hide()
        loss_chart.setTitle(self.tr("Loss vs Compression"))
        loss_chart.addSeries(loss_series)
        loss_chart.createDefaultAxes()
        loss_chart.axisX().setRange(0, 100)
        loss_chart.axisX().setTitleText(self.tr("quality (%)"))
        loss_chart.axisX().setTickCount(11)
        loss_chart.axisX().setLabelFormat("%d")
        loss_chart.axisY().setTitleText(self.tr("loss (%)"))
        loss_chart.setMinimumSize(600, 400)
        font = loss_chart.titleFont()
        font.setBold(True)
        loss_chart.setTitleFont(font)
        loss_view = QChartView(loss_chart)
        loss_view.setRenderHint(QPainter.Antialiasing)

        main_layout = QVBoxLayout()
        main_layout.addWidget(loss_view)
        self.setLayout(main_layout)
