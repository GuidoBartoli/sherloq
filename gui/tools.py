from PySide2.QtCore import Qt, Signal
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import (
    QTreeWidget,
    QHBoxLayout,
    QLabel,
    QTreeWidgetItem,
    QWidget, QSlider)

from utility import modify_font


class ToolWidget(QWidget):
    info_message = Signal(str)
    help_clicked = Signal(str)

    def __init__(self, parent=None):
        super(ToolWidget, self).__init__(parent)


class ToolTree(QTreeWidget):
    def __init__(self, parent=None):
        super(ToolTree, self).__init__(parent)
        group_names = []
        tool_names = []
        tool_infos = []
        tool_progress = []  # 0 = da fare, 1 = iniziato, 2 = funzionante, 3 = completo

        # [0]
        group_names.append(self.tr('[General]'))
        tool_names.append([self.tr('Original Image'),
                           self.tr('File Digest'),
                           self.tr('Hex Editor'),
                           self.tr('Reverse Search')])
        tool_infos.append([self.tr('Display the unaltered reference image for visual inspection'),
                           self.tr('Retrieve file information and compute many hashes and ballistics'),
                           self.tr('Open an external hexadecimal editor to show and edit raw bytes'),
                           self.tr('Use online search services to find visually similar images')])
        tool_progress.extend([3, 3, 2, 2])

        # [1]
        group_names.append(self.tr('[Metadata]'))
        tool_names.append([self.tr('Header Structure'),
                           self.tr('Metadata Extraction'),
                           self.tr('Thumbnail Analysis'),
                           self.tr('Geolocation data')])
        tool_infos.append([self.tr('Dump the physical EXIF structure and display an interactive view'),
                           self.tr('Scan through file metadata and gather all available information'),
                           self.tr('Extract optional embedded thumbnail and compare with original'),
                           self.tr('Retrieve optional geo-location data and show it on a world map')])
        tool_progress.extend([3, 3, 3, 2])

        # [2]
        group_names.append(self.tr('[Inspection]'))
        tool_names.append([self.tr('Enhancing Magnifier'),
                           self.tr('Reference Comparison'),
                           self.tr('Global Adjustments'),
                           self.tr('Fourier Transform')])
        tool_infos.append([self.tr('Use a loupe with visual enhancement for better identifying forgeries'),
                           self.tr('Open a synchronized double view to compare two different pictures'),
                           self.tr('Apply standard adjustments (contrast, brightness, hue, saturation)'),
                           self.tr('Compute amplitude and phase components of the 2D Fourier Transform')])
        tool_progress.extend([2, 2, 2, 3])

        # [3]
        group_names.append(self.tr('[JPEG]'))
        tool_names.append([self.tr('Quality Estimation'),
                           self.tr('Error Level Analysis'),
                           self.tr('Multiple Compression'),
                           self.tr('DCT Dimples Map')])
        tool_infos.append([self.tr('Extract quantization tables and estimate last saved JPEG quality'),
                           self.tr('Show pixel-level difference against different compression levels'),
                           self.tr('Use residuals to detect multiple compressions at different levels'),
                           self.tr('Analyze periodic quantization artifacts to detect manipulations')])
        tool_progress.extend([2, 3, 1, 0])

        # [4]
        group_names.append(self.tr('[Colors]'))
        tool_names.append([self.tr('RGB/HSV Plots'),
                           self.tr('PCA Projection'),
                           self.tr('Pixel Statistics'),
                           self.tr('Space Conversion')])
        tool_infos.append([self.tr('Display interactive 2D and 3D plots of RGB and HSV pixel data'),
                           self.tr('Use color PCA to project RGB values onto reduced vector spaces'),
                           self.tr('Compute minimum/maximum/average RGB values for every pixel'),
                           self.tr('Convert color channels into RGB/HSV/YCbCr/Lab/Luv/CMYK spaces')])
        tool_progress.extend([0, 2, 3, 2])

        # [5]
        group_names.append(self.tr('[Tonality]'))
        tool_names.append([self.tr('Luminance Gradient'),
                           self.tr('Echo Edge Filter'),
                           self.tr('Correlation Plot'),
                           self.tr('Wavelet Threshold')])
        tool_infos.append([self.tr('Analyze horizontal/vertical brightness variations across the image'),
                           self.tr('Use derivative filters to reveal artificial out-of-focus zones'),
                           self.tr('Exploit spatial correlation patterns among neighboring pixels'),
                           self.tr('Reconstruct image with different wavelet coefficient thresholds')])
        tool_progress.extend([2, 3, 0, 0])

        # [6]
        group_names.append(self.tr('[Noise]'))
        tool_names.append([self.tr('Noise Estimation'),
                           self.tr('Min/Max Deviation'),
                           self.tr('Image Bit Planes'),
                           self.tr('Frequency Separation')])
        tool_infos.append([self.tr('Estimate and visualize gaussian noise components of the image'),
                           self.tr('Highlight pixels deviating from block-based min/max statistics'),
                           self.tr('Visualize bit planes values to find different noise patterns'),
                           self.tr('Estimate high/low frequency components of the luminance channel')])
        tool_progress.extend([3, 2, 3, 0])

        # [7]
        group_names.append(self.tr('[Tampering]'))
        tool_names.append([self.tr('Contrast Enhancement'),
                           self.tr('Region Cloning'),
                           self.tr('Image Resampling'),
                           self.tr('Composite Splicing')])
        tool_infos.append([self.tr('Analyze color distribuions to detect contrast enhancements'),
                           self.tr('Use feature descriptors for copy/rotate clone area detection'),
                           self.tr('Analyze 2D pixel interpolation for detecting resampling traces'),
                           self.tr('Exploit DCT statistics for automatic splicing zone detection')])
        tool_progress.extend([0, 0, 0, 0])

        count = 0
        for i, group in enumerate(group_names):
            group_item = QTreeWidgetItem()
            group_item.setText(0, group)
            font = group_item.font(0)
            font.setBold(True)
            group_item.setFont(0, font)
            group_item.setData(0, Qt.UserRole, False)
            group_item.setIcon(0, QIcon('icons/{}.svg'.format(i)))
            for j, tool in enumerate(tool_names[i]):
                tool_item = QTreeWidgetItem(group_item)
                tool_item.setText(0, tool)
                tool_item.setData(0, Qt.UserRole, True)
                tool_item.setData(0, Qt.UserRole + 1, i)
                tool_item.setData(0, Qt.UserRole + 2, j)
                tool_item.setToolTip(0, tool_infos[i][j])
                if tool_progress[count] == 0:
                    modify_font(tool_item, italic=True)
                count += 1
            self.addTopLevelItem(group_item)
        self.expandAll()
        self.setColumnCount(1)
        self.header().setVisible(False)
        self.setMaximumWidth(300)
        self.version = '{:.2f}a'.format(sum(tool_progress) / (count * 3))

    def set_bold(self, tool, enabled):
        items = self.findItems(tool, Qt.MatchFixedString | Qt.MatchRecursive)
        if items:
            modify_font(items[0], bold=enabled)


class ParamSlider(QWidget):
    value_changed = Signal(int)

    def __init__(self, interval, step, ticks, reset=0, suffix='', parent=None):
        super(ParamSlider, self).__init__(parent)

        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(interval[0], interval[1])
        self.slider.setTickPosition(QSlider.TicksBelow)
        self.slider.setTickInterval((interval[1] - interval[0] + 1) / ticks)
        self.slider.setSingleStep(1)
        # self.slider.setPageStep(step)
        self.slider.setPageStep(1)
        self.slider.setValue(reset)
        self.slider.mouseDoubleClickEvent = self.double_click
        self.label = QLabel()
        modify_font(self.label, bold=True)
        self.suffix = suffix
        self.reset = reset
        self.sync(reset)
        self.slider.valueChanged.connect(self.sync)

        layout = QHBoxLayout()
        layout.addWidget(self.slider)
        layout.addWidget(self.label)
        self.setLayout(layout)

    def double_click(self, _):
        self.slider.setValue(self.reset)

    def sync(self, value):
        self.label.setText('{}{}'.format(value, self.suffix))
        self.value_changed.emit(value)

    def value(self):
        return self.slider.value()
