from PySide2.QtCharts import QtCharts
from PySide2.QtCore import Qt
from PySide2.QtGui import QPainter, QPen, QColor
from PySide2.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QCheckBox,
    QRadioButton,
    QWidget)


class HistWidget(QWidget):
    def __init__(self, parent=None):
        super(QWidget, self).__init__(parent)

        self.hist_chart = QtCharts.QChart()
        self.hist_chart.legend().hide()
        hist_view = QtCharts.QChartView(self.hist_chart)
        hist_view.setRenderHint(QPainter.Antialiasing)

        self.gray_radio = QRadioButton(self.tr('L'))
        self.red_radio = QRadioButton(self.tr('R'))
        self.green_radio = QRadioButton(self.tr('G'))
        self.blue_radio = QRadioButton(self.tr('B'))
        self.multi_radio = QRadioButton(self.tr('M'))
        self.log_check = QCheckBox(self.tr('Log'))

        self.gray_radio.toggled.connect(self.compute)
        self.red_radio.toggled.connect(self.compute)
        self.green_radio.toggled.connect(self.compute)
        self.blue_radio.toggled.connect(self.compute)
        self.multi_radio.toggled.connect(self.compute)
        self.log_check.stateChanged.connect(self.compute)

        button_layout = QHBoxLayout()
        button_layout.addWidget(self.gray_radio)
        button_layout.addWidget(self.red_radio)
        button_layout.addWidget(self.green_radio)
        button_layout.addWidget(self.blue_radio)
        button_layout.addWidget(self.multi_radio)
        button_layout.addWidget(self.log_check)
        button_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addWidget(hist_view)
        main_layout.addLayout(button_layout)
        self.setLayout(main_layout)
        self.image = None
        self.compute()

    def set_image(self, image):
        self.image = image
        self.compute()

    def compute(self):
        self.hist_chart.removeAllSeries()
        base_line = QtCharts.QLineSeries()
        for x in range(256):
            base_line.append(x, 0)
        self.hist_chart.addSeries(base_line)
        if self.image is not None:
            # if self.gray_radio.is
            # if self.gray:
            pass
            # hist_line = QtCharts.QLineSeries()
            # gray = cv.cvtColor(result, cv.COLOR_BGR2GRAY)
            # for i, h in enumerate([int(h[0]) for h in cv.calcHist([gray], [0], None, [256], [0, 256])]):
            #     hist_line.append(i, h)
            #     base_line.append(i, 0)
            # hist_line.setPen(QPen(Qt.white))
            # hist_area = QtCharts.QAreaSeries(hist_line, base_line)
            # color = QColor(Qt.lightGray)
            # color.setAlpha(128)
            # hist_area.setBrush(color)
            # hist_area.setPen(QPen(Qt.white))
            # self.hist_chart.addSeries(hist_line)

        # self.hist_chart.removeAllSeries()
        # hist_line = QtCharts.QLineSeries()
        # base_line = QtCharts.QLineSeries()
        # gray = cv.cvtColor(result, cv.COLOR_BGR2GRAY)
        # for i, h in enumerate([int(h[0]) for h in cv.calcHist([gray], [0], None, [256], [0, 256])]):
        #     hist_line.append(i, h)
        #     base_line.append(i, 0)
        # hist_line.setPen(QPen(Qt.white))
        # hist_area = QtCharts.QAreaSeries(hist_line, base_line)
        # color = QColor(Qt.lightGray)
        # color.setAlpha(128)
        # hist_area.setBrush(color)
        # hist_area.setPen(QPen(Qt.white))
        #
        # self.hist_chart.addSeries(hist_line)
        # self.hist_chart.createDefaultAxes()
        # self.hist_chart.axisX().setRange(0, 255)
        # self.hist_chart.axisX().setTitleText(self.tr('level'))
        # self.hist_chart.axisX().setTickCount(5)
        # self.hist_chart.axisX().setMinorTickCount(1)
        # self.hist_chart.axisX().setLabelFormat('%d')
        # # self.hist_chart.axisY().setRange(0, 100)
        # self.hist_chart.axisY().setTitleText(self.tr('count'))

        self.hist_chart.createDefaultAxes()
        self.hist_chart.axisX().setRange(0, 255)
        self.hist_chart.axisX().setTitleText(self.tr('level'))
        self.hist_chart.axisX().setTickCount(5)
        self.hist_chart.axisX().setMinorTickCount(1)
        self.hist_chart.axisX().setLabelFormat('%d')
        self.hist_chart.axisY().setRange(0, 100)
        self.hist_chart.axisY().setTitleText(self.tr('%'))




