from PySide6.QtCore import QUrl
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QVBoxLayout

from tools import ToolWidget


class EditorWidget(ToolWidget):
    def __init__(self, parent=None):
        super(EditorWidget, self).__init__(parent)
        web_view = QWebEngineView()
        web_view.load(QUrl("https://hexed.it/"))
        layout = QVBoxLayout()
        layout.addWidget(web_view)
        self.setLayout(layout)
