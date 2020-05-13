import cv2 as cv
import numpy as np
from PySide2.QtCore import Qt
from PySide2.QtWidgets import QVBoxLayout, QProgressDialog

from viewer import ImageViewer
from widget import ToolWidget


class GhostsWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(GhostsWidget, self).__init__(parent)


    def ghosts(self, image):
        gray = cv.cvtColor(image, cv.COLOR_BGR2GRAY).astype(float)

