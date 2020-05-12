import os
from subprocess import run, PIPE

from PySide2.QtCore import QUrl
from PySide2.QtWebEngineWidgets import QWebEngineView
from PySide2.QtWidgets import QVBoxLayout, QWidget


class StructureWidget(QWidget):
    def __init__(self, filename, parent=None):
        super(StructureWidget, self).__init__(parent)
        self.TEMP_FILE = os.path.join(os.getcwd(), 'structure.html')
        p = run(['pyexiftool/exiftool/exiftool', '-htmldump0', filename], stdout=PIPE)
        with open(self.TEMP_FILE, 'w') as file:
            file.write(p.stdout.decode('utf-8'))
        web_view = QWebEngineView()
        web_view.load(QUrl('file://' + self.TEMP_FILE))
        layout = QVBoxLayout()
        layout.addWidget(web_view)
        self.setLayout(layout)
        self.setMinimumWidth(750)

    def closeEvent(self, event):
        os.remove(self.TEMP_FILE)

