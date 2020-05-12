import cv2 as cv
from PySide2.QtCore import (
    QFileInfo,
    QLocale,
    QFile,
    QIODevice,
    QCryptographicHash)
from PySide2.QtWidgets import (
    QApplication,
    QAbstractItemView,
    QVBoxLayout,
    QWidget,
    QTableWidgetItem,
    QTableWidget,
    QMessageBox)

from utility import modify_font, human_size

