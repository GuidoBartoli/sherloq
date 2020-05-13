from PySide2.QtCore import Qt, Signal
from PySide2.QtGui import QIcon
from PySide2.QtWidgets import QTreeWidget, QTreeWidgetItem, QWidget

from utility import modify_font


class ToolWidget(QWidget):
    info_message = Signal(str)

    def __init__(self, parent=None):
        super(ToolWidget, self).__init__(parent)


class ToolTree(QTreeWidget):
    def __init__(self, parent=None):
        super(ToolTree, self).__init__(parent)
        group_names = []
        tool_names = []
        tool_infos = []

        group_names.append(self.tr('[General]'))                # [0]
        tool_names.append([self.tr('Original Image'),           # 3
                           self.tr('File Digest'),              # 2
                           self.tr('Similarity Search'),        # 0
                           self.tr('Automatic Tagging')])       # 0
        tool_infos.append([self.tr('Display the unaltered reference image for visual inspection'),
                           self.tr('Retrieve file information and compute many hashes and ballistics'),
                           self.tr('Use reverse search services for finding similar images online'),
                           self.tr('Apply deep learning algorithms for automatic picture tagging')])

        group_names.append(self.tr('[Metadata]'))               # [1]
        tool_names.append([self.tr('Header Structure'),         # 3
                           self.tr('Metadata Extraction'),      # 2
                           self.tr('Thumbnail Analysis'),       # 0
                           self.tr('Geolocation data')])        # 0
        tool_infos.append([self.tr('Dump the physical EXIF structure and display an interactive view'),
                           self.tr('Scan through file metadata and gather all available information'),
                           self.tr('Extract optional embedded thumbnail and compare with original'),
                           self.tr('Retrieve optional geo-location data and show it on a world map')])

        group_names.append(self.tr('[Inspection]'))             # [2]
        tool_names.append([self.tr('Enhancing Magnifier'),      # 0
                           self.tr('Reference Comparison'),     # 0
                           self.tr('Global Adjustments'),       # 0
                           self.tr('Range Compression')])       # 0
        tool_infos.append([self.tr('Use various visual enhancement for better identifying forgeries'),
                           self.tr('Open a synchronized double view to compare two different pictures'),
                           self.tr('Apply standard adjustments (contrast, brightness, hue, saturation)'),
                           self.tr('Compress tonality over different ranges to detect inconsistencies')])

        group_names.append(self.tr('[JPEG]'))                   # [3]
        tool_names.append([self.tr('Quality Estimation'),       # 0
                           self.tr('Error Level Analysis'),     # 3
                           self.tr('Quantization Ghosts'),      # 1
                           self.tr('Double Compression')])      # 0
        tool_infos.append([self.tr('Extract quantization tables and estimate last saved JPEG quality'),
                           self.tr('Show pixel-level difference against different compression levels'),
                           self.tr('Use residuals to detect multiple compressions at different levels'),
                           self.tr('Exploit DCT First-Digit-Statistics to detect double compression')])

        group_names.append(self.tr('[Colors]'))                 # [4]
        tool_names.append([self.tr('RGB/HSV 3D Plots'),         # 0
                           self.tr('RGB PCA Projection'),       # 0
                           self.tr('RGB Pixel Statistics'),     # 0
                           self.tr('Color Space Conversion')])  # 0
        tool_infos.append([self.tr('Display interactive 2D and 3D plots of RGB and HSV pixel data'),
                           self.tr('Use color PCA to project RGB values onto reduced vector spaces'),
                           self.tr('Compute Minimum/Maximum/Average RGB values for every pixel'),
                           self.tr('Convert color channels into RGB/HSV/YCbCr/Lab/CMYK color spaces')])

        group_names.append(self.tr('[Tonality]'))               # [5]
        tool_names.append([self.tr('Luminance Gradient'),       # 2
                           self.tr('Echo Edge Filter'),         # 0
                           self.tr('Frequency Separation'),     # 0
                           self.tr('Wavelet Reconstruction')])        # 0
        tool_infos.append([self.tr('Analyze horizontal and vertical brightness variations of the image'),
                           self.tr('Use derivative filter to reveal artificial out-of-focus zones'),
                           self.tr('Estimate high and low frequency components of the luminance channel'),
                           self.tr('Reconstruct image with different wavelet coefficient thresholds')])

        group_names.append(self.tr('[Noise]'))                  # [6]
        tool_names.append([self.tr('Noise Estimation'),         # 1
                           self.tr('Min/Max Deviation'),        # 0
                           self.tr('SNR Consistency'),          # 0
                           self.tr('Noise Segmentation')])      # 0
        tool_infos.append([self.tr('Estimate and visualize gaussian noise components of the image'),
                           self.tr('Highlight pixels deviating from block-based min/max statistics'),
                           self.tr('Evaluate uniformity of signal-to-noise ratio across the image'),
                           self.tr('Cluster noise into uniform regions for anomaly detection')])

        group_names.append(self.tr('[Tampering]'))              # [7]
        tool_names.append([self.tr('Contrast Enhancement'),     # 0
                           self.tr('Region Cloning'),           # 0
                           self.tr('Image Resampling'),         # 0
                           self.tr('Composite Splicing')])      # 0
        tool_infos.append([self.tr('Analyze color distribuions to detect contrast enhancements'),
                           self.tr('Use feature descriptors for copy/rotate clone area detection'),
                           self.tr('Analyze 2D pixel interpolation for detecting resampling traces'),
                           self.tr('Exploit DCT statistics for automatic splicing zone detection')])

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
            self.addTopLevelItem(group_item)
        self.expandAll()
        self.setColumnCount(1)
        self.header().setVisible(False)
        self.setMaximumWidth(300)

    def set_bold(self, tool_name, bold_enabled):
        items = self.findItems(tool_name, Qt.MatchFixedString | Qt.MatchRecursive)
        if items:
            modify_font(items[0], bold=bold_enabled)
