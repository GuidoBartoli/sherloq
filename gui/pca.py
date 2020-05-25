import cv2 as cv
import numpy as np
from PySide2.QtWidgets import (
    QAbstractItemView,
    QTableWidgetItem,
    QTableWidget,
    QVBoxLayout,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QRadioButton)

from tools import ToolWidget
from utility import normalize_mat, modify_font
from viewer import ImageViewer


class PcaWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(PcaWidget, self).__init__(parent)

        self.comp_combo = QComboBox()
        self.comp_combo.addItems([self.tr('#{}'.format(i + 1)) for i in range(3)])

        self.distvect_radio = QRadioButton(self.tr('Vector Distance'))
        self.cross_radio = QRadioButton(self.tr('Cross Correlation'))
        self.distvect_radio.setChecked(True)
        self.last_radio = self.distvect_radio

        self.image = image
        self.components = []
        rows, cols, dims = self.image.shape
        bgr = np.reshape(self.image, (rows * cols, dims)).astype(np.float32)
        m, eigen_vec, eigen_val = cv.PCACompute2(bgr, np.array([]))
        p = self.image.astype(np.float32) - m
        for v in eigen_vec:
            c = np.cross(p, v)
            d = np.linalg.norm(c, axis=2) / np.linalg.norm(v)
            distance = normalize_mat(d, to_bgr=True)
            cross = cv.merge([normalize_mat(x) for x in cv.split(c)])
            self.components.extend([distance, cross])

        table_data = [[m[0, 2], m[0, 1], m[0, 0]],
                      [eigen_vec[0, 2], eigen_vec[0, 1], eigen_vec[0, 0]],
                      [eigen_vec[1, 2], eigen_vec[1, 1], eigen_vec[1, 0]],
                      [eigen_vec[2, 2], eigen_vec[2, 1], eigen_vec[2, 0]],
                      [eigen_val[2, 0], eigen_val[1, 0], eigen_val[0, 0]]]
        table_widget = QTableWidget(5, 4)
        table_widget.setHorizontalHeaderLabels([
            self.tr('Element'), self.tr('Red'), self.tr('Green'), self.tr('Blue')])
        table_widget.setItem(0, 0, QTableWidgetItem(self.tr('Mean color')))
        table_widget.setItem(1, 0, QTableWidgetItem(self.tr('Eigen vect 1')))
        table_widget.setItem(2, 0, QTableWidgetItem(self.tr('Eigen vect 2')))
        table_widget.setItem(3, 0, QTableWidgetItem(self.tr('Eigen vect 3')))
        table_widget.setItem(4, 0, QTableWidgetItem(self.tr('Eigen values')))
        for i in range(len(table_data)):
            modify_font(table_widget.item(i, 0), bold=True)
            for j in range(len(table_data[i])):
                table_widget.setItem(i, j + 1, QTableWidgetItem(str(table_data[i][j])))
        # item = QTableWidgetItem()
        # item.setBackgroundColor(QColor(m[0, 2], m[0, 1], m[0, 0]))
        # table_widget.setItem(0, 4, item)
        # table_widget.resizeRowsToContents()
        # table_widget.resizeColumnsToContents()
        table_widget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table_widget.setSelectionMode(QAbstractItemView.SingleSelection)

        self.viewer = ImageViewer(self.image, self.image, None)
        self.process()

        self.comp_combo.currentIndexChanged.connect(self.process)
        self.distvect_radio.toggled.connect(self.process)
        self.cross_radio.toggled.connect(self.process)

        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr('Component:')))
        top_layout.addWidget(self.comp_combo)
        top_layout.addWidget(QLabel(self.tr('Projection:')))
        top_layout.addWidget(self.distvect_radio)
        top_layout.addWidget(self.cross_radio)
        top_layout.addStretch()
        bottom_layout = QHBoxLayout()
        bottom_layout.addWidget(table_widget)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        main_layout.addLayout(bottom_layout)
        self.setLayout(main_layout)

    def process(self):
        index = 2*self.comp_combo.currentIndex()
        if self.distvect_radio.isChecked():
            self.viewer.update_processed(self.components[index])
            self.last_radio = self.distvect_radio
        elif self.cross_radio.isChecked():
            self.viewer.update_processed(self.components[index + 1])
            self.last_radio = self.cross_radio
        else:
            self.last_radio.setChecked(True)
