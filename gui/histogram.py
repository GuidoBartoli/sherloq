import cv2 as cv
import numpy as np
from PySide2.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QCheckBox,
    QTableWidget,
    QTableWidgetItem,
    QAbstractItemView,
    QRadioButton)
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas)
from matplotlib.figure import Figure

from table import TableWidget
from tools import ToolWidget
from utility import compute_hist, modify_font


class HistWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(ToolWidget, self).__init__(parent)

        self.hist = [compute_hist(c) for c in cv.split(cv.cvtColor(image, cv.COLOR_BGR2RGB))]
        self.hist.append(compute_hist(cv.cvtColor(image, cv.COLOR_BGR2GRAY)))

        self.rgb_radio = QRadioButton(self.tr('RGB'))
        self.rgb_radio.setChecked(True)
        self.last_radio = self.rgb_radio
        self.red_radio = QRadioButton(self.tr('Red'))
        self.green_radio = QRadioButton(self.tr('Green'))
        self.blue_radio = QRadioButton(self.tr('Blue'))
        self.value_radio = QRadioButton(self.tr('Value'))
        self.log_check = QCheckBox(self.tr('Log scale'))
        self.grid_check = QCheckBox(self.tr('Show grid'))

        self.rgb_radio.toggled.connect(self.redraw)
        self.red_radio.toggled.connect(self.redraw)
        self.green_radio.toggled.connect(self.redraw)
        self.blue_radio.toggled.connect(self.redraw)
        self.value_radio.toggled.connect(self.redraw)
        self.log_check.stateChanged.connect(self.redraw)
        self.grid_check.stateChanged.connect(self.redraw)

        self.table_widget = QTableWidget(10, 2)
        self.table_widget.setHorizontalHeaderLabels([self.tr('Property'), self.tr('Value')])
        self.table_widget.setItem(0, 0, QTableWidgetItem(self.tr('Minimum')))
        self.table_widget.setItem(1, 0, QTableWidgetItem(self.tr('Average')))
        self.table_widget.setItem(2, 0, QTableWidgetItem(self.tr('Median')))
        self.table_widget.setItem(3, 0, QTableWidgetItem(self.tr('Maximum')))
        self.table_widget.setItem(4, 0, QTableWidgetItem(self.tr('Std dev')))
        self.table_widget.setItem(5, 0, QTableWidgetItem(self.tr('Pixels')))
        self.table_widget.setItem(6, 0, QTableWidgetItem(self.tr('Count')))
        self.table_widget.setItem(7, 0, QTableWidgetItem(self.tr('Percent')))
        self.table_widget.setItem(8, 0, QTableWidgetItem(self.tr('Unique')))
        self.table_widget.setItem(9, 0, QTableWidgetItem(self.tr('Smooth')))
        for i in range(self.table_widget.rowCount()):
            modify_font(self.table_widget.item(i, 0), bold=True)
        self.table_widget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table_widget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table_widget.resizeColumnsToContents()
        self.table_widget.setAlternatingRowColors(True)

        figure = Figure()
        plot_canvas = FigureCanvas(figure)
        self.axes = plot_canvas.figure.subplots()
        self.redraw()
        figure.set_tight_layout(True)

        center_layout = QHBoxLayout()
        center_layout.addWidget(plot_canvas)
        center_layout.addWidget(self.table_widget)

        bottom_layout = QHBoxLayout()
        bottom_layout.addWidget(self.rgb_radio)
        bottom_layout.addWidget(self.red_radio)
        bottom_layout.addWidget(self.green_radio)
        bottom_layout.addWidget(self.blue_radio)
        bottom_layout.addWidget(self.value_radio)
        bottom_layout.addWidget(self.log_check)
        bottom_layout.addWidget(self.grid_check)
        bottom_layout.addStretch()

        main_layout = QVBoxLayout()
        main_layout.addLayout(center_layout)
        main_layout.addLayout(bottom_layout)
        self.setLayout(main_layout)

    def redraw(self):
        # FIXME: la funzione viene chiamata due volte ad ogni evento
        x = np.arange(256)
        alpha = 0.25
        self.axes.clear()
        if self.value_radio.isChecked():
            y = self.hist[3]
            self.axes.step(x, y, 'k')
            self.axes.fill_between(x, y, alpha=alpha, facecolor='k', step='pre')
        else:
            if self.red_radio.isChecked() or self.rgb_radio.isChecked():
                y = self.hist[0]
                self.axes.step(x, y, 'r')
                self.axes.fill_between(x, y, alpha=alpha, facecolor='r', step='pre')
            if self.green_radio.isChecked() or self.rgb_radio.isChecked():
                y = self.hist[1]
                self.axes.step(x, y, 'g')
                self.axes.fill_between(x, y, alpha=alpha, facecolor='g', step='pre')
            if self.blue_radio.isChecked() or self.rgb_radio.isChecked():
                y = self.hist[2]
                self.axes.step(x, y, 'b')
                self.axes.fill_between(x, y, alpha=alpha, facecolor='b', step='pre')
        if self.log_check.isChecked():
            self.axes.set_yscale('log')
            self.axes.set_ylim(bottom=1)
        else:
            self.axes.set_yscale('linear')
            self.axes.set_ylim(bottom=0)
        self.axes.set_xlim([-1, 256])
        self.axes.grid(self.grid_check.isChecked(), which='both')
        self.axes.figure.canvas.draw()

        if self.rgb_radio.isChecked():
            self.table_widget.setEnabled(False)
        else:
            _, _, (_, argmin), (_, argmax) = cv.minMaxLoc(y)
            pixels = np.sum(y)
            mean = np.sum(x * y) / pixels
            stddev = np.sqrt(np.sum(((x - mean)**2) * y) / pixels)
            median = 0
            self.table_widget.setItem(0, 1, QTableWidgetItem(str(argmin)))
            self.table_widget.setItem(1, 1, QTableWidgetItem(str(np.round(mean, 2))))
            self.table_widget.setItem(2, 1, QTableWidgetItem(str(median)))
            self.table_widget.setItem(3, 1, QTableWidgetItem(str(argmax)))
            self.table_widget.setItem(4, 1, QTableWidgetItem(str(np.round(stddev, 2))))
            self.table_widget.setItem(5, 1, QTableWidgetItem(str(pixels)))
            self.table_widget.resizeColumnsToContents()
            self.table_widget.setEnabled(True)
