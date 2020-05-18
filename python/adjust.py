import cv2 as cv
import numpy as np
from PySide2.QtCore import Qt
from PySide2.QtGui import QPainter, QPen, QColor
from PySide2.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QGridLayout,
    QCheckBox,
    QPushButton,
    QComboBox)
from PySide2.QtCharts import QtCharts

from slider import ParamSlider
from tools import ToolWidget
from utility import create_lut, signed_value, modify_font
from viewer import ImageViewer


class AdjustWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(AdjustWidget, self).__init__(parent)

        params_layout = QGridLayout()
        params_layout.addWidget(QLabel(self.tr('Brightness:')), 0, 0)
        self.bright_slider = ParamSlider([-255, +255], [1, 8], 16)
        self.bright_slider.valueChanged.connect(self.process)
        params_layout.addWidget(self.bright_slider, 0, 1)
        self.bright_label = QLabel()
        modify_font(self.bright_label, bold=True)
        params_layout.addWidget(self.bright_label, 0, 2)

        params_layout.addWidget(QLabel(self.tr('Saturation:')), 1, 0)
        self.sat_slider = ParamSlider([-255, +255], [1, 8], 16)
        self.sat_slider.valueChanged.connect(self.process)
        params_layout.addWidget(self.sat_slider, 1, 1)
        self.sat_label = QLabel()
        modify_font(self.sat_label, bold=True)
        params_layout.addWidget(self.sat_label, 1, 2)

        params_layout.addWidget(QLabel(self.tr('Hue:')), 2, 0)
        self.hue_slider = ParamSlider([0, 180], [1, 5], 10)
        self.hue_slider.valueChanged.connect(self.process)
        params_layout.addWidget(self.hue_slider, 2, 1)
        self.hue_label = QLabel()
        modify_font(self.hue_label, bold=True)
        params_layout.addWidget(self.hue_label, 2, 2)

        params_layout.addWidget(QLabel(self.tr('Gamma:')), 0, 3)
        self.gamma_slider = ParamSlider([0, 50], [1, 5], 10)
        self.gamma_slider.setValue(10)
        self.gamma_slider.valueChanged.connect(self.process)
        params_layout.addWidget(self.gamma_slider, 0, 4)
        self.gamma_label = QLabel()
        modify_font(self.gamma_label, bold=True)
        params_layout.addWidget(self.gamma_label, 0, 5)

        params_layout.addWidget(QLabel(self.tr('Shadows:')), 1, 3)
        self.shadow_slider = ParamSlider([-100, +100], [1, 10], 10)
        self.shadow_slider.valueChanged.connect(self.process)
        params_layout.addWidget(self.shadow_slider, 1, 4)
        self.shadow_label = QLabel()
        modify_font(self.shadow_label, bold=True)
        params_layout.addWidget(self.shadow_label, 1, 5)

        params_layout.addWidget(QLabel(self.tr('Highlights:')), 2, 3)
        self.high_slider = ParamSlider([-100, +100], [1, 10], 10)
        self.high_slider.valueChanged.connect(self.process)
        params_layout.addWidget(self.high_slider, 2, 4)
        self.high_label = QLabel()
        modify_font(self.high_label, bold=True)
        params_layout.addWidget(self.high_label, 2, 5)

        self.equalize_combo = QComboBox()
        self.equalize_combo.addItems(
            [self.tr('No equalization'), self.tr('Histogram EQ'), self.tr('Weak CLAHE'),
             self.tr('Medium CLAHE'), self.tr('Strong CLAHE'), self.tr('Extreme CLAHE')])
        self.equalize_combo.currentIndexChanged.connect(self.process)
        params_layout.addWidget(self.equalize_combo, 0, 6)
        self.invert_check = QCheckBox(self.tr('Invert'))
        self.invert_check.stateChanged.connect(self.process)
        params_layout.addWidget(self.invert_check, 1, 6)
        self.reset_button = QPushButton(self.tr('Reset'))
        self.reset_button.clicked.connect(self.reset)
        params_layout.addWidget(self.reset_button, 2, 6)

        self.hist_chart = QtCharts.QChart()
        self.hist_chart.setPlotAreaBackgroundBrush(Qt.lightGray)
        self.hist_chart.setPlotAreaBackgroundVisible(True)
        self.hist_chart.setBackgroundVisible(False)
        self.hist_chart.legend().hide()
        hist_view = QtCharts.QChartView(self.hist_chart)
        hist_view.setRenderHint(QPainter.Antialiasing)
        top_layout = QHBoxLayout()
        top_layout.addLayout(params_layout)
        top_layout.addWidget(hist_view)

        self.image = image
        self.viewer = ImageViewer(self.image, self.image)
        self.process()

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def process(self):
        brightness = self.bright_slider.value()
        self.bright_label.setText(signed_value(brightness))
        saturation = self.sat_slider.value()
        self.sat_label.setText(signed_value(saturation))
        hue = self.hue_slider.value()
        self.hue_label.setText('{}°'.format(2*hue))
        gamma = self.gamma_slider.value() / 10
        self.gamma_label.setText('{:.1f}'.format(gamma))
        shadows = self.shadow_slider.value()
        self.shadow_label.setText('{}%'.format(signed_value(shadows)))
        highlights = self.high_slider.value()
        self.high_label.setText('{}%'.format(signed_value(highlights)))
        equalize = self.equalize_combo.currentIndex()
        invert = self.invert_check.isChecked()

        result = np.copy(self.image)
        if brightness != 0 or saturation != 0 or hue != 0:
            h, s, v = cv.split(cv.cvtColor(result, cv.COLOR_BGR2HSV))
            if hue != 0:
                h = h.astype(np.float64) + hue
                h[h < 0] += 180
                h[h > 180] -= 180
                h = h.astype(np.uint8)
            if saturation != 0:
                s = cv.add(s, saturation)
            if brightness != 0:
                v = cv.add(v, brightness)
            result = cv.cvtColor(cv.merge((h, s, v)), cv.COLOR_HSV2BGR)
        if gamma != 0:
            inverse = 1 / gamma
            lut = np.array([((i / 255) ** inverse) * 255 for i in np.arange(0, 256)]).astype(np.uint8)
            result = cv.LUT(result, lut)
        if shadows != 0:
            result = cv.LUT(result, create_lut(int(shadows / 100 * 255), 0))
        if highlights != 0:
            result = cv.LUT(result, create_lut(0, int(highlights / 100 * 255)))
        if equalize > 0:
            h, s, v = cv.split(cv.cvtColor(result, cv.COLOR_BGR2HSV))
            if equalize == 1:
                v = cv.equalizeHist(v)
            elif equalize > 1:
                clip = 0
                if equalize == 2:
                    clip = 2
                elif equalize == 3:
                    clip = 10
                elif equalize == 4:
                    clip = 20
                elif equalize == 5:
                    clip = 40
                v = cv.createCLAHE(clip).apply(v)
            result = cv.cvtColor(cv.merge((h, s, v)), cv.COLOR_HSV2BGR)
        if invert:
            result = cv.bitwise_not(result)
        self.viewer.update_processed(result)

        self.hist_chart.removeAllSeries()
        hist_line = QtCharts.QLineSeries()
        base_line = QtCharts.QLineSeries()
        gray = cv.cvtColor(result, cv.COLOR_BGR2GRAY)
        for i, h in enumerate([int(h[0]) for h in cv.calcHist([gray], [0], None, [256], [0, 256])]):
            hist_line.append(i, h)
            base_line.append(i, 0)
        hist_line.setPen(QPen(Qt.white))
        hist_area = QtCharts.QAreaSeries(hist_line, base_line)
        color = QColor(Qt.lightGray)
        color.setAlpha(128)
        hist_area.setBrush(color)
        hist_area.setPen(QPen(Qt.white))

        self.hist_chart.addSeries(hist_line)
        self.hist_chart.createDefaultAxes()
        self.hist_chart.axisX().setRange(0, 255)
        self.hist_chart.axisX().setTitleText(self.tr('level'))
        self.hist_chart.axisX().setTickCount(5)
        self.hist_chart.axisX().setMinorTickCount(1)
        self.hist_chart.axisX().setLabelFormat('%d')
        # self.hist_chart.axisY().setRange(0, 100)
        self.hist_chart.axisY().setTitleText(self.tr('count'))

    def reset(self):
        self.bright_slider.setValue(0)
        self.sat_slider.setValue(0)
        self.hue_slider.setValue(0)
        self.gamma_slider.setValue(10)
        self.shadow_slider.setValue(0)
        self.high_slider.setValue(0)
        self.equalize_combo.setCurrentIndex(0)
        self.invert_check.setChecked(False)
