from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem, QWidget

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
        tool_progress = []  # 0 = da fare, 1 = debug, 2 = funzionante, 3 = completo

        # [0]
        group_names.append(self.tr("[General]"))
        tool_names.append(
            [
                self.tr("Original Image"),
                self.tr("File Digest"),
                self.tr("Hex Editor"),
                self.tr("Similar Search"),
            ]
        )
        tool_infos.append(
            [
                self.tr("Display the unaltered reference image for visual inspection"),
                self.tr(
                    "Retrieve physical file information, crypto and perceptual hashes"
                ),
                self.tr(
                    "Open an external hexadecimal editor to show and edit raw bytes"
                ),
                self.tr(
                    "Browse online search services to find visually similar images"
                ),
            ]
        )
        tool_progress.extend([3, 3, 2, 2])

        # [1]
        group_names.append(self.tr("[Metadata]"))
        tool_names.append(
            [
                self.tr("Header Structure"),
                self.tr("EXIF Full Dump"),
                self.tr("Thumbnail Analysis"),
                self.tr("Geolocation data"),
            ]
        )
        tool_infos.append(
            [
                self.tr(
                    "Dump the file header structure and display an interactive view"
                ),
                self.tr(
                    "Scan through file metadata and gather all available information"
                ),
                self.tr(
                    "Extract optional embedded thumbnail and compare with original"
                ),
                self.tr(
                    "Retrieve optional geolocation data and show it on a world map"
                ),
            ]
        )
        tool_progress.extend([3, 3, 3, 2])

        # [2]
        group_names.append(self.tr("[Inspection]"))
        tool_names.append(
            [
                self.tr("Enhancing Magnifier"),
                self.tr("Channel Histogram"),
                self.tr("Global Adjustments"),
                self.tr("Reference Comparison"),
            ]
        )
        tool_infos.append(
            [
                self.tr(
                    "Magnifying glass with enhancements for better identifying forgeries"
                ),
                self.tr(
                    "Display single color channels or RGB composite interactive histogram"
                ),
                self.tr(
                    "Apply standard image adjustments (brightness, hue, saturation, ...)"
                ),
                self.tr(
                    "Open a synchronized double view for comparison with another picture"
                ),
            ]
        )
        tool_progress.extend([3, 3, 3, 3])

        # [3]
        group_names.append(self.tr("[Detail]"))
        tool_names.append(
            [
                self.tr("Luminance Gradient"),
                self.tr("Echo Edge Filter"),
                self.tr("Wavelet Threshold"),
                self.tr("Frequency Split"),
            ]
        )
        tool_infos.append(
            [
                self.tr(
                    "Analyze horizontal/vertical brightness variations across the image"
                ),
                self.tr(
                    "Use derivative filters to reveal artificial out-of-focus regions"
                ),
                self.tr(
                    "Reconstruct image with different wavelet coefficient thresholds"
                ),
                self.tr(
                    "Divide image luminance into high and low frequency components"
                ),
            ]
        )
        tool_progress.extend([3, 3, 3, 3])

        # [4]
        group_names.append(self.tr("[Colors]"))
        tool_names.append(
            [
                self.tr("RGB/HSV Plots"),
                self.tr("Space Conversion"),
                self.tr("PCA Projection"),
                self.tr("Pixel Statistics"),
            ]
        )
        tool_infos.append(
            [
                self.tr(
                    "Display interactive 2D and 3D plots of RGB and HSV pixel values"
                ),
                self.tr("Convert RGB channels into HSV/YCbCr/Lab/Luv/CMYK/Gray spaces"),
                self.tr("Use color PCA to project pixel onto most salient components"),
                self.tr("Compute minimum/maximum/average RGB values for every pixel"),
            ]
        )
        tool_progress.extend([3, 3, 3, 3])

        # [5]
        group_names.append(self.tr("[Noise]"))
        tool_names.append(
            [
                self.tr("Signal Separation"),
                self.tr("Min/Max Deviation"),
                self.tr("Bit Plane Values"),
                self.tr("Wavelet Blocking"),
                self.tr("PRNU Identification"),
            ]
        )
        tool_infos.append(
            [
                self.tr(
                    "Estimate and extract different kind of image noise components"
                ),
                self.tr(
                    "Highlight pixels deviating from block-based min/max statistics"
                ),
                self.tr(
                    "Show individual bit planes to find inconsistent noise patterns"
                ),
                self.tr("Noise estimation based on high pass wavelet coefficients & grid blocking"),
                self.tr("Exploit sensor pattern noise introduced by different cameras"),
            ]
        )
        tool_progress.extend([3, 3, 3, 2, 0])

        # [6]
        group_names.append(self.tr("[JPEG]"))
        tool_names.append(
            [
                self.tr("Quality Estimation"),
                self.tr("Error Level Analysis"),
                self.tr("Multiple Compression"),
                self.tr("JPEG Ghost Maps"),
            ]
        )
        tool_infos.append(
            [
                self.tr(
                    "Extract quantization tables and estimate last saved JPEG quality"
                ),
                self.tr(
                    "Show pixel-wise differences against a fixed compression level"
                ),
                self.tr("Use a machine learning model to detect multiple compression"),
                self.tr(
                    "Highlight traces of different compressions in difference images"
                ),
            ]
        )
        tool_progress.extend([3, 3, 0, 2])

        # [7]
        group_names.append(self.tr("[Tampering]"))
        tool_names.append(
            [
                self.tr("Contrast Enhancement"),
                self.tr("Copy-Move Forgery"),
                self.tr("Composite Splicing"),
                self.tr("Image Resampling"),
            ]
        )
        tool_infos.append(
            [
                self.tr("Analyze color distributions to detect contrast enhancements"),
                self.tr("Use invariant feature descriptors to detect cloned regions"),
                self.tr("Exploit DCT statistics for automatic splicing zone detection"),
                self.tr(
                    "Estimate 2D pixel interpolation for detecting resampling traces"
                ),
            ]
        )
        tool_progress.extend([3, 2, 3, 2])


        # [8]
        group_names.append(self.tr("[AI Solutions]"))
        tool_names.append(
            [
                self.tr("TruFor"),
            ]
        )
        tool_infos.append(
            [
                self.tr("TruFor: Leveraging all-round clues for trustworthy image forgery detection and localization"),
                
            ]
        )
        tool_progress.extend([2])


        # [9]
        group_names.append(self.tr("[Various]"))
        tool_names.append(
            [
                self.tr("Median Filtering"),
                self.tr("Illuminant Map"),
                self.tr("Dead/Hot Pixels"),
                self.tr("Stereogram Decoder"),
            ]
        )
        tool_infos.append(
            [
                self.tr("Detect nonlinear processing traces left by median filtering"),
                self.tr(
                    "Estimate scene local light direction on estimated 3D surfaces"
                ),
                self.tr(
                    "Detect and fix dead/hot pixels caused by sensor imperfections"
                ),
                self.tr(
                    "Decode 3D images concealed inside crossed-eye autostereograms"
                ),
            ]
        )
        tool_progress.extend([2, 0, 0, 3])

        count = 0
        for i, group in enumerate(group_names):
            group_item = QTreeWidgetItem()
            group_item.setText(0, group)
            font = group_item.font(0)
            font.setBold(True)
            group_item.setFont(0, font)
            group_item.setData(0, Qt.UserRole, False)
            group_item.setIcon(0, QIcon(f"icons/{i}.svg"))
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
        self.version = f"{sum(tool_progress) / 100:.2f}f"

    def set_bold(self, tool, enabled):
        items = self.findItems(tool, Qt.MatchFixedString | Qt.MatchRecursive)
        if items:
            modify_font(items[0], bold=enabled)
