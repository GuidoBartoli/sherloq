import subprocess
import subprocess
from shutil import copyfile

import numpy as np
from PySide2.QtCore import QTemporaryFile
from PySide2.QtWidgets import (
    QSizePolicy,
    QPlainTextEdit,
    QLabel,
    QGridLayout,
    QMessageBox)
from tabulate import tabulate

from jpeg import TABLE_SIZE, ZIG_ZAG, DCT_SIZE, get_tables
from tools import ToolWidget
from utility import modify_font, get_exiftool
from table import TableWidget


class QualityWidget(ToolWidget):
    def __init__(self, filename, parent=None):
        super(QualityWidget, self).__init__(parent)

        MRK = b'\xFF'
        SOI = b'\xD8'
        DQT = b'\xDB'
        MSK = b'\x0F'
        PAD = b'\x00'

        MAX_TABLES = 2
        LEN_OFFSET = 2
        LUMA_IDX = 0
        CHROMA_IDX = 1

        temp_file = QTemporaryFile()
        if temp_file.open():
            copyfile(filename, temp_file.fileName())
            subprocess.run([get_exiftool(), '-all=', '-overwrite_original', temp_file.fileName()],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            luma = np.zeros((DCT_SIZE, DCT_SIZE), dtype=int)
            chroma = np.zeros((DCT_SIZE, DCT_SIZE), dtype=int)
            found = False
            with open(temp_file.fileName(), 'rb') as file:
                first = file.read(1)
                if first not in [MRK, SOI]:
                    QMessageBox.warning(self, self.tr('Warning'), self.tr('File is not a JPEG image!'))
                    return
                while True:
                    if not self.find_next(file, [MRK, DQT, PAD]):
                        break
                    length = file.read(1)[0] - LEN_OFFSET
                    if length <= 0 or length % (TABLE_SIZE + 1) != 0:
                        continue
                    while length > 0:
                        mode = file.read(1)
                        if not mode:
                            break
                        index = mode[0] & MSK[0]
                        if index >= MAX_TABLES:
                            break
                        length -= 1
                        for k in range(TABLE_SIZE):
                            b = file.read(1)[0]
                            if not b:
                                break
                            length -= 1
                            i, j = ZIG_ZAG[k]
                            if index == LUMA_IDX:
                                luma[i, j] = b
                            elif index == CHROMA_IDX:
                                chroma[i, j] = b
                        else:
                            found = True
        if not found:
            QMessageBox.warning(self, self.tr('Warning'), self.tr('Unable to find JPEG tables!'))
            return

        tables = np.concatenate((luma[:, :, np.newaxis], chroma[:, :, np.newaxis]), axis=2)
        levels = [0, 0]
        for i in range(2):
            table = tables[:, :, i]
            mean = np.mean(table) * TABLE_SIZE
            mean = (mean - table[0, 0]) / (TABLE_SIZE - 1)
            levels[i] = (1 - mean / 255) * 100
        profile = [np.mean(np.abs(tables - get_tables(q))) for q in range(100)]
        quality = np.argmin(profile)

        # headers = [str(i) for i in range(DCT_SIZE)]
        # luma_widget = TableWidget(luma.tolist(), headers, bold=False, align=True, search=False)

        luma_text = QPlainTextEdit()
        modify_font(luma_text, mono=True)
        luma_text.setReadOnly(True)
        luma_text.appendPlainText(self.tr('Luminance QT (level = {:.2f}%)\n'.format(levels[0])))
        luma_text.appendPlainText(tabulate(luma, tablefmt='plain'))

        chroma_text = QPlainTextEdit()
        chroma_text.setReadOnly(True)
        modify_font(chroma_text, mono=True)
        chroma_text.appendPlainText(self.tr('Chrominance QT (level = {:.2f}%)\n'.format(levels[1])))
        chroma_text.appendPlainText(tabulate(chroma, tablefmt='plain'))

        quality_label = QLabel(self.tr('Estimated JPEG quality (last save) = {}%'.format(quality)))
        modify_font(quality_label, bold=True)
        deviation_label = QLabel(self.tr('(deviation from standard tables = {:.2f})'.format(profile[quality])))
        modify_font(deviation_label, italic=True)

        main_layout = QGridLayout()
        main_layout.addWidget(luma_text, 0, 0)
        main_layout.addWidget(chroma_text, 0, 1)
        main_layout.addWidget(quality_label, 1, 0)
        main_layout.addWidget(deviation_label, 1, 1)
        self.setLayout(main_layout)
        self.setMinimumSize(590, 270)

    def find_next(self, file, markers):
        while True:
            for m in markers:
                b = file.read(1)
                if not b:
                    return False
                if b != m:
                    break
            else:
                return True
