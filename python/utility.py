from PySide2.QtGui import QImage, QFontDatabase
import numpy as np
from PySide2.QtWidgets import QTreeWidgetItem
import cv2 as cv


def mat2img(cvmat):
    height, width, channel = cvmat.shape
    return QImage(cvmat.data, width, height, 3 * width, QImage.Format_BGR888)


def modify_font(obj, bold=False, italic=False, underline=False, mono=False):
    if mono:
        font = QFontDatabase.systemFont(QFontDatabase.FixedFont)
    else:
        if type(obj) is QTreeWidgetItem:
            font = obj.font(0)
        else:
            font = obj.font()
    font.setBold(bold)
    font.setItalic(italic)
    font.setUnderline(underline)
    if type(obj) is QTreeWidgetItem:
        obj.setFont(0, font)
    else:
        obj.setFont(font)


def human_size(total, binary=False, suffix='B'):
    units = ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']
    if binary:
        units = [unit + 'i' for unit in units]
        factor = 1024.0
    else:
        factor = 1000.0
    for unit in units:
        if abs(total) < factor:
            return "%3.1f %s%s" % (total, unit, suffix)
        total /= factor
    return "%.1f %s%s" % (total, units[-1], suffix)


def compress_jpeg(image, quality):
    result, buffer = cv.imencode(".jpg", image, [cv.IMWRITE_JPEG_QUALITY, quality])
    if not result:
        raise RuntimeError('Utility.compress_jpg: Unable to encode JPG!')
    return cv.imdecode(buffer, cv.IMREAD_COLOR)


def create_lut(low, high):
    if low >= 0:
        p1 = (+low, 0)
    else:
        p1 = (0, -low)
    if high >= 0:
        p2 = (255 - high, 255)
    else:
        p2 = (255, 255 + high)
    lut = [(x * (p1[1] - p2[1]) + p1[0] * p2[1] - p1[1] * p2[0]) / (p1[0] - p2[0]) for x in range(256)]
    return np.clip(np.array(lut), 0, 255).astype(np.uint8)
