import cv2 as cv
from PySide2.QtWidgets import (
    QVBoxLayout)

from viewer import ImageViewer
from widget import ToolWidget


class NoiseWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super(NoiseWidget, self).__init__(parent)
        self.image = image
        self.viewer = ImageViewer(image, image)
        self.process()
        layout = QVBoxLayout()
        layout.addWidget(self.viewer)
        self.setLayout(layout)

    def process(self):
        gray = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        denoised = cv.medianBlur(gray, 3)
        noise = cv.normalize(cv.absdiff(gray.astype(float), denoised.astype(float)), None, 0, 255)
        self.viewer.update_processed(cv.cvtColor(noise.astype(int), cv.COLOR_GRAY2BGR))
