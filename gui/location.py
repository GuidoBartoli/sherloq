import os

from PySide2.QtCore import QUrl, QTemporaryDir, Qt
from PySide2.QtWebEngineWidgets import QWebEngineView
from PySide2.QtWidgets import QVBoxLayout, QLabel

from pyexiftool import exiftool
from tools import ToolWidget
from utility import exiftool_exe, modify_font


class LocationWidget(ToolWidget):
    def __init__(self, filename, parent=None):
        super(LocationWidget, self).__init__(parent)
        self.temp_dir = QTemporaryDir()
        if self.temp_dir.isValid():
            with exiftool.ExifTool(exiftool_exe()) as et:
                try:
                    metadata = et.get_metadata(filename)
                    lat = metadata["Composite:GPSLatitude"]
                    lon = metadata["Composite:GPSLongitude"]
                except KeyError:
                    label = QLabel(self.tr("Geolocation data not found!"))
                    modify_font(label, bold=True)
                    label.setStyleSheet("color: #FF0000")
                    label.setAlignment(Qt.AlignCenter)
                    layout = QVBoxLayout()
                    layout.addWidget(label)
                    self.setLayout(layout)
                    return
                url = (
                    f"https://www.google.com/maps/place/{lat},{lon}/@{lat},{lon},17z/"
                    f"data=!4m5!3m4!1s0x0:0x0!8m2!3d{lat}!4d{lon}"
                )
                web_view = QWebEngineView()
                web_view.load(QUrl(url))
                layout = QVBoxLayout()
                layout.addWidget(web_view)
                self.setLayout(layout)
