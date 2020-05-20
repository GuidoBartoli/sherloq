from PySide2.QtCore import Qt, Signal
from PySide2.QtWidgets import QSlider


class ParamSlider(QSlider):
    def __init__(self, interval, steps, ticks, reset, parent=None):
        super(ParamSlider, self).__init__(Qt.Horizontal, parent)

        self.setRange(interval[0], interval[1])
        self.setTickPosition(QSlider.TicksBelow)
        self.setTickInterval((interval[1] - interval[0] + 1) / ticks)
        self.setSingleStep(steps[0])
        self.setPageStep(steps[1])
        self.reset = reset

    def mouseDoubleClickEvent(self, event):
        self.setValue(self.reset)
