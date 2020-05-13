import os

import cv2 as cv
import magic
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
    QTableWidgetItem,
    QTableWidget,
    QMessageBox)

from utility import modify_font, human_size
from widget import ToolWidget


class DigestWidget(ToolWidget):
    def __init__(self, filename, image, parent=None):
        super(DigestWidget, self).__init__(parent)
        rows = 17
        cols = 2

        table_widget = QTableWidget()
        table_widget.setColumnCount(cols)
        table_widget.setRowCount(rows)
        headers = [self.tr('Property'), self.tr('Value')]
        table_widget.setHorizontalHeaderLabels(headers)
        modify_font(table_widget, mono=True)

        file_info = QFileInfo(filename)
        table_widget.setItem(0, 0, QTableWidgetItem(self.tr('File name')))
        table_widget.setItem(0, 1, QTableWidgetItem(file_info.fileName()))
        table_widget.setItem(1, 0, QTableWidgetItem(self.tr('Parent folder')))
        table_widget.setItem(1, 1, QTableWidgetItem(str(file_info.dir().absolutePath())))
        table_widget.setItem(2, 0, QTableWidgetItem(self.tr('File type')))
        table_widget.setItem(2, 1, QTableWidgetItem(magic.from_file(filename, mime=True)))
        table_widget.setItem(3, 0, QTableWidgetItem(self.tr('File size')))
        table_widget.setItem(3, 1, QTableWidgetItem('{} bytes ({})'.format(
            QLocale().toString(file_info.size()), human_size(file_info.size()))))
        table_widget.setItem(4, 0, QTableWidgetItem(self.tr('File owner')))
        table_widget.setItem(4, 1, QTableWidgetItem(file_info.owner()))
        table_widget.setItem(5, 0, QTableWidgetItem(self.tr('Permissions')))
        table_widget.setItem(5, 1, QTableWidgetItem(str(oct(os.stat(filename).st_mode)[-3:])))
        table_widget.setItem(6, 0, QTableWidgetItem(self.tr('Creation time')))
        table_widget.setItem(6, 1, QTableWidgetItem(file_info.birthTime().toLocalTime().toString()))
        table_widget.setItem(7, 0, QTableWidgetItem(self.tr('Last access')))
        table_widget.setItem(7, 1, QTableWidgetItem(file_info.lastRead().toLocalTime().toString()))
        table_widget.setItem(8, 0, QTableWidgetItem(self.tr('Last modified')))
        table_widget.setItem(8, 1, QTableWidgetItem(file_info.lastModified().toLocalTime().toString()))
        table_widget.setItem(9, 0, QTableWidgetItem(self.tr('Metadata change')))
        table_widget.setItem(9, 1, QTableWidgetItem(file_info.metadataChangeTime().toLocalTime().toString()))

        # width = image.shape[1]
        # height = image.shape[0]
        # pixels = human_size(width*height, suffix='px')
        # table_widget.setItem(8, 0, QTableWidgetItem(self.tr('Resolution')))
        # table_widget.setItem(8, 1, QTableWidgetItem('{}x{} px ({})'.format(height, width, pixels)))
        # table_widget.setItem(6, 0, QTableWidgetItem(self.tr('Channels')))
        # table_widget.setItem(6, 1, QTableWidgetItem('{} (RGB)'.format(image.shape[2])))
        # progress = QProgressDialog(self.tr('Counting unique colors...'), None, 0, width*height, self)
        # progress.setWindowModality(Qt.WindowModal)
        # rgb_hist = np.zeros((256, 256, 256), int)
        # i = 0
        # for r in range(height):
        #     for c in range(width):
        #         blue, green, red = image[r][c]
        #         rgb_hist[red][green][blue] += 1
        #         progress.setValue(i)
        #         i += 1
        # table_widget.setItem(4, 0, QTableWidgetItem(self.tr('Colors')))
        # table_widget.setItem(4, 1, QTableWidgetItem(str(np.count_nonzero(rgb_hist))))

        file = QFile(filename)
        if not file.open(QIODevice.ReadOnly):
            QMessageBox.warning(self, self.tr('Warning'), self.tr('Unable to read file from disk!'))
            return
        data = file.readAll()
        md5 = QCryptographicHash.hash(data, QCryptographicHash.Md5).toHex()
        table_widget.setItem(10, 0, QTableWidgetItem(self.tr('MD5')))
        table_widget.setItem(10, 1, QTableWidgetItem(str(md5, encoding='utf-8')))
        sha1 = QCryptographicHash.hash(data, QCryptographicHash.Sha1).toHex()
        table_widget.setItem(11, 0, QTableWidgetItem(self.tr('SHA2-1')))
        table_widget.setItem(11, 1, QTableWidgetItem(str(sha1, encoding='utf-8')))
        sha224 = QCryptographicHash.hash(data, QCryptographicHash.Sha224).toHex()
        table_widget.setItem(12, 0, QTableWidgetItem(self.tr('SHA2-224')))
        table_widget.setItem(12, 1, QTableWidgetItem(str(sha224, encoding='utf-8')))
        sha256 = QCryptographicHash.hash(data, QCryptographicHash.Sha256).toHex()
        table_widget.setItem(13, 0, QTableWidgetItem(self.tr('SHA2-256')))
        table_widget.setItem(13, 1, QTableWidgetItem(str(sha256, encoding='utf-8')))
        sha3_256 = QCryptographicHash.hash(data, QCryptographicHash.Sha3_256).toHex()
        table_widget.setItem(14, 0, QTableWidgetItem(self.tr('SHA3-256')))
        table_widget.setItem(14, 1, QTableWidgetItem(str(sha3_256, encoding='utf-8')))

        phash = cv.img_hash.pHash(image)
        table_widget.setItem(15, 0, QTableWidgetItem(self.tr('pHash')))
        table_widget.setItem(15, 1, QTableWidgetItem(str(phash[0])))
        avg_hash = cv.img_hash.averageHash(image)
        table_widget.setItem(16, 0, QTableWidgetItem(self.tr('aHash')))
        table_widget.setItem(16, 1, QTableWidgetItem(str(avg_hash[0])))

        for row in range(rows):
            modify_font(table_widget.item(row, 0), bold=True)
        table_widget.resizeColumnsToContents()
        table_widget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        layout = QVBoxLayout()
        layout.addWidget(table_widget)
        self.setLayout(layout)
        self.adjustSize()
        self.setMinimumSize(700, 570)
        table_widget.itemDoubleClicked.connect(self.copy_cell)

    def copy_cell(self, item):
        QApplication.clipboard().setText(item.text())
        self.message_to_show.emit(self.tr('Cell contents copied to clipboard'))

