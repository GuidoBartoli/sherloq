import cv2 as cv
from time import time
import numpy as np
from PySide2.QtWidgets import (
    QSpinBox,
    QCheckBox,
    QComboBox,
    QLabel,
    QVBoxLayout,
    QDoubleSpinBox,
    QGridLayout)
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure

from tools import ToolWidget
from utility import elapsed_time


class PlotsWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(PlotsWidget, self).__init__(parent)

        choices = ['Red', 'Green', 'Blue', 'Hue', 'Saturation', 'Value']
        self.xaxis_combo = QComboBox()
        self.xaxis_combo.addItems(choices)
        self.xaxis_combo.setCurrentIndex(3)
        self.yaxis_combo = QComboBox()
        self.yaxis_combo.addItems(choices)
        self.yaxis_combo.setCurrentIndex(4)
        self.sampling_spin = QSpinBox()
        levels = int(np.log2(min(image.shape[:-1])))
        self.sampling_spin.setRange(0, levels)
        self.sampling_spin.setSpecialValueText(self.tr('Off'))
        self.sampling_spin.setValue(1)
        self.sampling_spin.setSuffix(self.tr(' level(s)'))
        self.size_spin = QSpinBox()
        self.size_spin.setRange(1, 10)
        self.size_spin.setValue(2)
        self.size_spin.setSuffix(self.tr(' pt'))
        self.style_combo = QComboBox()
        self.markers = [',', '.', 'o', '8', 's', 'p', 'P', '*', 'h', 'H', 'X', 'D']
        self.style_combo.addItems(
            ['pixel', 'point', 'circle', 'octa', 'square', 'penta', 'plus', 'star', 'hexa1', 'hexa2', 'cross', 'diamond'])
        self.alpha_spin = QDoubleSpinBox()
        self.alpha_spin.setRange(0, 1)
        self.alpha_spin.setDecimals(2)
        self.alpha_spin.setSingleStep(0.05)
        self.alpha_spin.setValue(1)
        self.colors_check = QCheckBox(self.tr('Show colors'))
        self.grid_check = QCheckBox(self.tr('Show grid'))
        self.norm_check = QCheckBox(self.tr('Normalized'))
        self.total_label = QLabel()

        img = np.copy(image)
        self.colors = [None] * (levels + 1)
        for scale in range(levels + 1):
            rgb = cv.cvtColor(img.astype(np.float32) / 255, cv.COLOR_BGR2RGB)
            hsv = cv.cvtColor(rgb, cv.COLOR_RGB2HSV)
            hsv[:, :, 0] /= 360
            shape = (img.shape[0] * img.shape[1], img.shape[2])
            self.colors[scale] = np.concatenate((np.reshape(rgb, shape), np.reshape(hsv, shape)), axis=1)
            img = cv.pyrDown(img)
        figure = Figure()
        plot_canvas = FigureCanvas(figure)
        self.axes = plot_canvas.figure.subplots()
        self.redraw()
        figure.set_tight_layout(True)

        self.xaxis_combo.currentIndexChanged.connect(self.redraw)
        self.yaxis_combo.currentIndexChanged.connect(self.redraw)
        self.sampling_spin.valueChanged.connect(self.redraw)
        self.size_spin.valueChanged.connect(self.redraw)
        self.style_combo.currentIndexChanged.connect(self.redraw)
        self.alpha_spin.valueChanged.connect(self.redraw)
        self.colors_check.stateChanged.connect(self.redraw)
        self.grid_check.stateChanged.connect(self.redraw)
        self.norm_check.stateChanged.connect(self.redraw)

        bottom_layout = QGridLayout()
        bottom_layout.addWidget(NavigationToolbar(plot_canvas, self), 0, 0, 1, 5)
        bottom_layout.addWidget(QLabel(self.tr('X axis:')), 1, 0)
        bottom_layout.addWidget(self.xaxis_combo, 1, 1)
        bottom_layout.addWidget(QLabel(self.tr('Y Axis:')), 2, 0)
        bottom_layout.addWidget(self.yaxis_combo, 2, 1)
        bottom_layout.addWidget(QLabel(self.tr('Subsampling:')), 3, 0)
        bottom_layout.addWidget(self.sampling_spin, 3, 1)
        bottom_layout.addWidget(QLabel(self.tr('Point size:')), 1, 2)
        bottom_layout.addWidget(self.size_spin, 1, 3)
        bottom_layout.addWidget(QLabel(self.tr('Point style:')), 2, 2)
        bottom_layout.addWidget(self.style_combo, 2, 3)
        bottom_layout.addWidget(QLabel(self.tr('Point alpha:')), 3, 2)
        bottom_layout.addWidget(self.alpha_spin, 3, 3)
        bottom_layout.addWidget(self.colors_check, 1, 4)
        bottom_layout.addWidget(self.grid_check, 2, 4)
        bottom_layout.addWidget(self.total_label, 3, 4)

        main_layout = QVBoxLayout()
        main_layout.addWidget(plot_canvas)
        main_layout.addLayout(bottom_layout)
        self.setLayout(main_layout)

    def redraw(self):
        start = time()
        v = self.sampling_spin.value()
        x = self.colors[v][:, self.xaxis_combo.currentIndex()]
        y = self.colors[v][:, self.yaxis_combo.currentIndex()]
        xlim = self.axes.get_xlim()
        ylim = self.axes.get_ylim()
        s = self.size_spin.value()**2
        c = None if not self.colors_check.isChecked() else self.colors[v][:, :3]
        a = self.alpha_spin.value()
        m = self.markers[self.style_combo.currentIndex()]
        self.axes.clear()
        self.axes.set_facecolor([0.5] * 3 if c is not None else [1.0] * 3)
        self.axes.scatter(x, y, s, c, m, alpha=a)
        self.axes.set_xlabel(self.xaxis_combo.currentText())
        self.axes.set_ylabel(self.yaxis_combo.currentText())
        self.axes.grid(self.grid_check.isChecked(), which='both')
        self.axes.set_xlim(xlim)
        self.axes.set_ylim(ylim)
        self.axes.figure.canvas.draw()
        self.total_label.setText(self.tr('[{} points]'.format(len(x))))
        self.info_message.emit('Plot redraw = {}'.format(elapsed_time(start)))
