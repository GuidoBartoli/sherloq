import os
import sys
from subprocess import run, PIPE

from PySide6.QtCore import QUrl, QTemporaryDir
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QVBoxLayout

from tools import ToolWidget
from utility import exiftool_exe


class HeaderWidget(ToolWidget):
    def __init__(self, filename, parent=None):
        super(HeaderWidget, self).__init__(parent)
        self.temp_dir = QTemporaryDir()
        if self.temp_dir.isValid():
            temp_file = os.path.join(self.temp_dir.path(), "structure.html")
            input_data = b"\n"
            command = [exiftool_exe(), "-htmldump0", filename]
            p = run(command, input=input_data, stdout=PIPE)
            with open(temp_file, "w") as file:
                file.write(p.stdout.decode("utf-8"))

            web_view = QWebEngineView()
            if sys.platform.startswith("win32"):
                temp_file = self.temp_dir.path() + "//" + "structure.html"  # fixes the broken ESC sequence
                web_view.load(QUrl(temp_file))  # load without 'file://' prefix
            else:
                web_view.load(QUrl("file://" + temp_file))

            layout = QVBoxLayout()
            layout.addWidget(web_view)
            self.setLayout(layout)
            self.setMinimumWidth(900)
