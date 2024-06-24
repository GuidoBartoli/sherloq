import subprocess
from shutil import copyfile

import cv2 as cv
import numpy as np
from PySide6.QtCore import QTemporaryFile, Qt
from PySide6.QtGui import QColor, QBrush
from joblib import load
from PySide6.QtWidgets import (
    QLabel,
    QVBoxLayout,
    QGridLayout,
    QTableWidget,
    QMessageBox,
    QTableWidgetItem,
    QAbstractItemView,
)
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

from jpeg import TABLE_SIZE, ZIG_ZAG, DCT_SIZE, get_tables, loss_curve, compress_jpg
from tools import ToolWidget
from utility import modify_font, exiftool_exe, clip_value


class QualityWidget(ToolWidget):
    def __init__(self, filename, image, parent=None):
        super(QualityWidget, self).__init__(parent)

        x = np.arange(1, 101)
        y = loss_curve(image)
        tail = 5
        qm = np.argmin(y[:-tail]) + 1
        if qm == 100 - tail:
            qm = 100

        figure = Figure()
        canvas = FigureCanvas(figure)
        axes = canvas.figure.subplots()
        axes.plot(x, y * 100, label="compression loss")
        axes.fill_between(x, y * 100, alpha=0.2)
        axes.axvline(qm, linestyle=":", color="k", label=f"min error (q = {qm})")
        xt = axes.get_xticks()
        xt = np.append(xt, 1)
        axes.set_xticks(xt)
        axes.set_xlim([1, 100])
        axes.set_ylim([0, 100])
        axes.set_xlabel(self.tr("JPEG quality (%)"))
        axes.set_ylabel(self.tr("average error (%)"))
        axes.grid(True, which="both")
        axes.legend(loc="upper center")
        axes.figure.canvas.draw()
        figure.set_tight_layout(True)

        main_layout = QVBoxLayout()
        main_layout.addWidget(canvas)

        MRK = b"\xFF"
        SOI = b"\xD8"
        DQT = b"\xDB"
        # DHT = b'\xC4'
        MSK = b"\x0F"
        PAD = b"\x00"

        MAX_TABLES = 2
        LEN_OFFSET = 2
        LUMA_IDX = 0
        CHROMA_IDX = 1

        luma = np.zeros((DCT_SIZE, DCT_SIZE), dtype=int)
        chroma = np.zeros((DCT_SIZE, DCT_SIZE), dtype=int)
        temp_file = QTemporaryFile()
        try:
            if temp_file.open():
                copyfile(filename, temp_file.fileName())
                subprocess.run(
                    [
                        exiftool_exe(),
                        "-all=",
                        "-overwrite_original",
                        temp_file.fileName(),
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                found = False
                with open(temp_file.fileName(), "rb") as file:
                    first = file.read(1)
                    if first not in [MRK, SOI]:
                        raise ValueError(self.tr("File is not a JPEG image!"))
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
                raise ValueError(self.tr("Unable to find JPEG tables!"))

            levels = [
                (1 - (np.mean(t.ravel()[1:]) - 1) / 254) * 100 for t in [luma, chroma]
            ]
            distance = np.zeros(101)
            for qm in range(101):
                lu, ch = cv.split(get_tables(qm))
                lu_diff = np.mean(cv.absdiff(luma, lu))
                ch_diff = np.mean(cv.absdiff(chroma, ch))
                distance[qm] = (lu_diff + 2 * ch_diff) / 3
            closest = np.argmin(distance)
            deviation = distance[closest]
            if deviation == 0:
                quality = closest
                message = "(standard tables)"
            else:
                quality = int(np.round(closest - deviation))
                message = f"(deviation from standard tables --> {deviation:.4f})"
            if quality == 0:
                quality = 1
            quality_label = QLabel(
                self.tr(f"[JPEG FORMAT] Last saved quality: {quality}% {message}")
            )
            modify_font(quality_label, bold=True)

            luma_label = QLabel(
                self.tr(f"Luminance Quantization Table (level = {levels[0]:.2f}%)\n")
            )
            luma_label.setAlignment(Qt.AlignCenter)
            modify_font(luma_label, underline=True)
            luma_table = self.create_table(luma)
            luma_table.setFixedSize(420, 190)

            chroma_label = QLabel(
                self.tr(f"Chrominance Quantization Table (level = {levels[1]:.2f}%)\n")
            )
            chroma_label.setAlignment(Qt.AlignCenter)
            modify_font(chroma_label, underline=True)
            chroma_table = self.create_table(chroma)
            chroma_table.setFixedSize(420, 190)

            table_layout = QGridLayout()
            table_layout.addWidget(luma_label, 0, 0)
            table_layout.addWidget(luma_table, 1, 0)
            table_layout.addWidget(chroma_label, 0, 1)
            table_layout.addWidget(chroma_table, 1, 1)
            table_layout.addWidget(quality_label, 2, 0, 1, 2)
            main_layout.addLayout(table_layout)

        except ValueError:
            modelfile = "models/jpeg_qf.mdl"
            try:
                model = load(modelfile)
                limit = (
                    model.best_ntree_limit
                    if hasattr(model, "best_ntree_limit")
                    else None
                )
                # f = self.get_features(image)
                # p = model.predict_proba(f, ntree_limit=limit)[0, 0]
                qp = model.predict(np.reshape(y, (1, len(y))), ntree_limit=limit)[0]
                # if p > 0.5:
                #     p = 2 * (p - 0.5) * 100
                #     output = self.tr('Uncompressed image (p = {:.2f}%)'.format(p))
                # else:
                #     p = (1 - 2 * p) * 100
                #     output = self.tr('Compressed image (p = {:.2f}%) ---> Estimated JPEG quality = {}%'.format(p, qm))
                message = self.tr(
                    f"[LOSSLESS FORMAT] Estimated last saved quality = {qp:.1f}%{'' if qp <= 99 else ' (uncompressed)'}"
                )
                if qp == 100:
                    message += " (uncompressed)"
                prob_label = QLabel(message)
                modify_font(prob_label, bold=True)
                main_layout.addWidget(prob_label)
            except FileNotFoundError:
                QMessageBox.critical(
                    self, self.tr("Error"), self.tr(f'Model not found ("{modelfile}")!')
                )

        main_layout.addStretch()
        self.setLayout(main_layout)

    def show_error(self, message):
        error_label = QLabel(message)
        modify_font(error_label, bold=True)
        error_label.setStyleSheet("color: #FF0000")
        error_label.setAlignment(Qt.AlignCenter)
        main_layout = QVBoxLayout()
        main_layout.addWidget(error_label)
        self.setLayout(main_layout)

    @staticmethod
    def get_features(image):
        # q = [6, 8, 10, 13, 31, 71, 73, 75, 76, 81, 84, 87, 88, 93, 94, 96, 97, 98, 99, 100]
        q = list(range(1, 101))
        c = loss_curve(image, q)
        return np.reshape(c, (1, len(q)))

    @staticmethod
    def find_next(file, markers):
        while True:
            for m in markers:
                b = file.read(1)
                if not b:
                    return False
                if b != m:
                    break
            else:
                return True

    @staticmethod
    def create_table(matrix):
        table_widget = QTableWidget(DCT_SIZE, DCT_SIZE)
        hsv = np.array([[[0, 192, 255]]])
        maximum = clip_value(np.max(matrix) - 1, minv=1)
        for i in range(DCT_SIZE):
            for j in range(DCT_SIZE):
                value = matrix[i, j]
                item = QTableWidgetItem(str(value))
                item.setTextAlignment(Qt.AlignCenter)
                hsv[0, 0, 0] = 64 - 64 * ((value - 1) / maximum)
                rgb = cv.cvtColor(hsv.astype(np.uint8), cv.COLOR_HSV2RGB)
                item.setBackground(
                    QBrush(QColor(rgb[0, 0, 0], rgb[0, 0, 1], rgb[0, 0, 2]))
                )
                table_widget.setItem(i, j, item)
        table_widget.resizeRowsToContents()
        table_widget.resizeColumnsToContents()
        table_widget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        table_widget.setSelectionMode(QAbstractItemView.SingleSelection)
        modify_font(table_widget, mono=True)
        return table_widget
