from PySide2.QtCore import Signal
from PySide2.QtWidgets import QWidget


class ToolWidget(QWidget):
    message_to_show = Signal(str)

    def __init__(self, parent=None):
        super(ToolWidget, self).__init__(parent)
