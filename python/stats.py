import cv2 as cv
import numpy as np
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QCheckBox,
    QLabel,
    QRadioButton,
    QPushButton)

from tools import ToolWidget
from viewer import ImageViewer


class StatsWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(StatsWidget, self).__init__(parent)

        params_layout = QHBoxLayout()
        params_layout.addWidget(QLabel(self.tr('Mode:')))
        self.min_radio = QRadioButton(self.tr('Minimum'))
        self.min_radio.setChecked(True)
        self.min_radio.toggled.connect(self.process)
        params_layout.addWidget(self.min_radio)
        self.avg_radio = QRadioButton(self.tr('Average'))
        self.avg_radio.toggled.connect(self.process)
        params_layout.addWidget(self.avg_radio)
        self.max_radio = QRadioButton(self.tr('Maximum'))
        self.max_radio.toggled.connect(self.process)
        params_layout.addWidget(self.max_radio)
        self.incl_check = QCheckBox(self.tr('Inclusive'))
        self.incl_check.stateChanged.connect(self.process)
        params_layout.addWidget(self.incl_check)
        params_layout.addStretch()
        help_button = QPushButton()
        help_button.setIcon(QIcon('icons/help.svg'))
        help_button.clicked.connect(self.show_help)
        params_layout.addWidget(help_button)

        self.image = image
        self.viewer = ImageViewer(self.image, self.image)
        self.process()

        main_layout = QVBoxLayout()
        main_layout.addLayout(params_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def process(self):
        inclusive = self.incl_check.isChecked()
        b, g, r = cv.split(self.image)
        if self.min_radio.isChecked():
            if inclusive:
                blue_mask = np.logical_and(b <= g, b <= r)
                green_mask = np.logical_and(g <= r, g <= b)
                red_mask = np.logical_and(r <= b, r <= g)
            else:
                blue_mask = np.logical_and(b < g, b < r)
                green_mask = np.logical_and(g < r, g < b)
                red_mask = np.logical_and(r < b, r < g)
        elif self.max_radio.isChecked():
            if inclusive:
                blue_mask = np.logical_and(b >= g, b >= r)
                green_mask = np.logical_and(g >= r, g >= b)
                red_mask = np.logical_and(r >= b, r >= g)
            else:
                blue_mask = np.logical_and(b > g, b > r)
                green_mask = np.logical_and(g > r, g > b)
                red_mask = np.logical_and(r > b, r > g)
        elif self.avg_radio.isChecked():
            if inclusive:
                blue_mask = np.logical_or(np.logical_and(b >= r, b <= g), np.logical_and(b >= g, b <= r))
                green_mask = np.logical_or(np.logical_and(g >= r, g <= b), np.logical_and(g >= b, g <= r))
                red_mask = np.logical_or(np.logical_and(r >= b, r <= g), np.logical_and(r >= g, r <= b))
            else:
                blue_mask = np.logical_or(np.logical_and(b > r, b < g), np.logical_and(b > g, b < r))
                green_mask = np.logical_or(np.logical_and(g > r, g < b), np.logical_and(g > b, g < r))
                red_mask = np.logical_or(np.logical_and(r > b, r < g), np.logical_and(r > g, r < b))
        else:
            return
        result = np.zeros_like(self.image)
        result[blue_mask] = np.array([255, 0, 0])
        result[green_mask] = np.array([0, 255, 0])
        result[red_mask] = np.array([0, 0, 255])
        self.viewer.update_processed(result)

    def show_help(self):
        self.help_clicked.emit('RGB Pixel Statistics')
