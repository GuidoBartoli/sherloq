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

        group_names.append(self.tr('[General]'))
        tool_names.append([self.tr('Original Image'),
                           self.tr('File Digest'),
                           self.tr('Similarity Search'),
                           self.tr('Automatic Tagging')])
        tool_infos.append([self.tr('Display the unaltered reference image for visual inspection'),
                           self.tr('Compute various hashes together with extension ballistics'),
                           self.tr('Use reverse search services for finding similar images online'),
                           self.tr('Apply deep learning algorithms for automatic picture tagging')])

        group_names.append(self.tr('[Metadata]'))
        tool_names.append([self.tr('Header Structure'),
                           self.tr('Metadata Extraction'),
                           self.tr('Thumbnail Analysis'),
                           self.tr('Geolocation data')])
        tool_infos.append([self.tr('Dump the physical EXIF structure and display an interactive view'),
                           self.tr('Scan through file metadata and gather all available information'),
                           self.tr('Extract embedded thumbnail (if present) and compare with original'),
                           self.tr('Retrieve geo-location data (if present) and show on a world map')])

        group_names.append(self.tr('[Inspection]'))
        tool_names.append([self.tr('Enhancing Magnifier'),
                           self.tr('Reference Comparison'),
                           self.tr('Global Adjustments'),
                           self.tr('Range Compression')])
        tool_infos.append([self.tr('Use various visual enhancement for better identifying forgeries'),
                           self.tr('Open a synchronized double view to compare two different pictures'),
                           self.tr('Apply standard adjustments (contrast, brightness, hue, saturation)'),
                           self.tr('Compress tonality over different ranges to detect inconsistencies')])

        group_names.append(self.tr('[JPEG]'))
        tool_names.append([self.tr('Quality Estimation'),
                           self.tr('Error Level Analysis'),
                           self.tr('Quantization Ghosts'),
                           self.tr('Double Compression')])
        tool_infos.append([self.tr('Extract quantization tables and estimate last saved JPEG quality'),
                           self.tr('Show pixel-level difference against different compression levels'),
                           self.tr('Use residuals to detect multiple compressions at different levels'),
                           self.tr('Exploit DCT First-Digit-Statistics to detect double compression')])

        group_names.append(self.tr('[Colors]'))
        tool_names.append([self.tr('RGB/HSV 3D Plots'),
                           self.tr('RGB PCA Projection'),
                           self.tr('RGB Pixel Statistics'),
                           self.tr('Color Space Conversion')])
        tool_infos.append([self.tr('Display interactive 2D/3D plots of RGB and HSV pixel data'),
                           self.tr('Use color PCA to project RGB values onto a different vector space'),
                           self.tr('Compute Minimum/Maximum/Average RGB values for every pixel'),
                           self.tr('Convert color channels into RGB/HSV/YCbCr/Lab/CMYK color spaces')])

        group_names.append(self.tr('[Tonality]'))
        tool_names.append([self.tr('Luminance Gradient'),
                           self.tr('Echo Edge Filter'),
                           self.tr('Frequency Separation'),
                           self.tr('Wavelet Analysis')])
        tool_infos.append([self.tr('Analyze horizontal and vertical brightness variations of the image'),
                           self.tr('Estimate the high-frequency component of the luminance channel'),
                           self.tr('Use derivative filter to reveal artificial out-of-focus zones'),
                           self.tr('Reconstruct image with different wavelet coefficient thresholds')])

        group_names.append(self.tr('[Noise]'))
        tool_names.append([self.tr('Noise Estimation'),
                           self.tr('Min/Max Deviation'),
                           self.tr('SNR Consistency'),
                           self.tr('Noise Segmentation')])
        tool_infos.append([self.tr('Estimate and visualize gaussian noise components of the image'),
                           self.tr('Highlight pixels deviating from block-based min/max statistics'),
                           self.tr('Evaluate uniformity of signal-to-noise ratio across the image'),
                           self.tr('Cluster noise into uniform regions for anomaly detection')])

        group_names.append(self.tr('[Tampering]'))
        tool_names.append([self.tr('Contrast Enhancement'),
                           self.tr('Region Cloning'),
                           self.tr('Image Resampling'),
                           self.tr('Composite Splicing')])
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
