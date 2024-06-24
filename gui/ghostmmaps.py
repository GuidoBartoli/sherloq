# This code implements JPEG Ghost maps as explained in the paper: "Exposing Digital Forgeries from JPEG Ghosts" by Hany Farid
# The book "Digital Image Forensics" by Hany Farid gives a more detailed explanation of the technique for those interested

from PySide6.QtWidgets import (
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QSpinBox,
    QCheckBox,
    QPushButton,
)

from tools import ToolWidget
from viewer import ImageViewer

# ghost necessary imports
import matplotlib.pyplot as plt
import math
import numpy as np
import cv2
import os


class GhostmapWidget(ToolWidget):
    # tool layout
    def __init__(self, filename, image, parent=None):
        super(GhostmapWidget, self).__init__(parent)

        # save variables to self
        self.filename = filename

        # store different xy-offsets so user can quickly cycle different maps and inspect changes
        self.ghostmaps = [None] * 64

        # prepare user interface - input variables
        # qmin
        self.qmin_spin = QSpinBox()
        self.qmin_spin.setRange(0, 100)
        self.qmin_spin.setValue(50)
        # qmax
        self.qmax_spin = QSpinBox()
        self.qmax_spin.setRange(0, 100)
        self.qmax_spin.setValue(90)
        # qstep
        self.qstep_spin = QSpinBox()
        self.qstep_spin.setRange(0, 20)
        self.qstep_spin.setValue(5)
        # lattice X offset
        self.xoffset_spin = QSpinBox()
        self.xoffset_spin.setRange(0, 7)
        self.xoffset_spin.setValue(0)
        # lattice Y offset
        self.yoffset_spin = QSpinBox()
        self.yoffset_spin.setRange(0, 7)
        self.yoffset_spin.setValue(0)
        # grayscale
        self.showgray_check = QCheckBox(self.tr("Grayscale"))
        self.showgray_check.setChecked(True)
        # plot original above the maps?
        self.includeoriginal_check = QCheckBox(self.tr("Include original in plot?"))
        self.includeoriginal_check.setChecked(False)
        # calculate ghost maps
        self.process_button = QPushButton(self.tr("Calculate Maps"))
        # calculate next offset ghost maps
        self.process_next_offset_button = QPushButton(self.tr("Next offset"))
        # calculate previous offset ghost maps
        self.process_previous_offset_button = QPushButton(self.tr("Previous offset"))

        # combine top layout
        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr("Lower Quality:")))
        top_layout.addWidget(self.qmin_spin)
        top_layout.addWidget(QLabel(self.tr("Upper Quality:")))
        top_layout.addWidget(self.qmax_spin)
        top_layout.addWidget(QLabel(self.tr("Quality Step:")))
        top_layout.addWidget(self.qstep_spin)
        top_layout.addWidget(self.showgray_check)
        top_layout.addWidget(self.includeoriginal_check)
        top_layout.addWidget(self.process_button)
        top_layout.addWidget(QLabel(self.tr("Offset X:")))
        top_layout.addWidget(self.xoffset_spin)
        top_layout.addWidget(QLabel(self.tr("Offset Y:")))
        top_layout.addWidget(self.yoffset_spin)
        top_layout.addWidget(self.process_previous_offset_button)
        top_layout.addWidget(self.process_next_offset_button)
        top_layout.addStretch()

        self.viewer = ImageViewer(image, image, None)

        self.plt = plt
        self.processGhostmaps()

        self.process_button.clicked.connect(self.processGhostmaps)
        self.process_previous_offset_button.clicked.connect(
            self.calculate_previous_offset
        )
        self.process_next_offset_button.clicked.connect(self.calculate_next_offset)

        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)

    def calculate_next_offset(self):
        x_offset = self.xoffset_spin.value()
        y_offset = self.yoffset_spin.value()

        if x_offset < 7:
            x_offset = x_offset + 1
            self.xoffset_spin.setValue(x_offset)
            self.processGhostmaps()
        elif y_offset < 7:
            x_offset = 0
            y_offset = y_offset + 1
            self.xoffset_spin.setValue(x_offset)
            self.yoffset_spin.setValue(y_offset)
            self.processGhostmaps()
        else:
            x_offset = 0
            y_offset = 0
            self.xoffset_spin.setValue(x_offset)
            self.yoffset_spin.setValue(y_offset)
            self.processGhostmaps()

    def calculate_previous_offset(self):
        x_offset = self.xoffset_spin.value()
        y_offset = self.yoffset_spin.value()

        if x_offset > 0:
            x_offset = x_offset - 1
            self.xoffset_spin.setValue(x_offset)
            self.processGhostmaps()
        elif y_offset > 0:
            x_offset = 7
            y_offset = y_offset - 1
            self.xoffset_spin.setValue(x_offset)
            self.yoffset_spin.setValue(y_offset)
            self.processGhostmaps()
        else:
            x_offset = 7
            y_offset = 7
            self.xoffset_spin.setValue(x_offset)
            self.yoffset_spin.setValue(y_offset)
            self.processGhostmaps()

    # calculate ghost maps fuction:
    def processGhostmaps(self):
        self.process_button.setEnabled(False)  # wait for processing
        self.plt.clf()

        Qmin = self.qmin_spin.value()
        Qmax = self.qmax_spin.value()
        Qstep = self.qstep_spin.value()

        averagingBlock = 16

        shift_x = self.xoffset_spin.value()
        shift_y = self.yoffset_spin.value()

        includeoriginal = self.includeoriginal_check.isChecked()
        grayscale = self.showgray_check.isChecked()

        # if ghostmaps already exists, return map and exit
        if self.ghostmaps[shift_x + shift_y * 8] is not None:
            (
                ghostplot,
                sqmin,
                sqmax,
                sqstep,
                sincludeoriginal,
                sgrayscale,
            ) = self.ghostmaps[shift_x + shift_y * 8]
            if (
                sqmin == Qmin
                and sqmax == Qmax
                and sqstep == Qstep
                and sincludeoriginal == includeoriginal
                and sgrayscale == grayscale
            ):
                self.viewer.update_processed(ghostplot)
                self.process_button.setEnabled(True)
                return

        # load original
        original = np.double(cv2.imread(self.filename))
        ydim, xdim, zdim = original.shape

        # construct ghostmaps with possible shift
        nQ = int((Qmax - Qmin) / Qstep) + 1

        i = 0
        ghostmap = np.zeros((ydim, xdim, nQ))
        for quality in range(Qmin, Qmax + 1, Qstep):
            # Shift the image because:
            # misalignment of JPEG block lattice may destroy the JPEG ghost since new spatial frequencies will be introduced
            # by shifting we can search for the correct allignment, if there is one
            shifted_original = np.roll(original, shift_x, axis=1)
            shifted_original = np.roll(shifted_original, shift_y, axis=0)
            # compute difference original and re-compressed versions of original
            # store recompressed image as variable so we're not writing to disk:
            tempvar1 = cv2.imencode(
                ".jpg", shifted_original, [int(cv2.IMWRITE_JPEG_QUALITY), quality]
            )[1].tobytes()
            tempcar2 = np.frombuffer(tempvar1, np.byte)
            tmpResave = np.double(cv2.imdecode(tempcar2, cv2.IMREAD_ANYCOLOR))

            # compute difference and average over RGB
            for z in range(zdim):
                ghostmap[:, :, i] += np.square(
                    shifted_original[:, :, z].astype(np.double) - tmpResave[:, :, z]
                )

            ghostmap[:, :, i] /= zdim
            i += 1

        # compute average over larger area to counter complicating factor, as explained in paper
        blkE = np.zeros(
            (int((ydim) / averagingBlock), int((xdim) / averagingBlock), nQ)
        )
        for c in range(nQ):
            cy = 0
            for y in range(0, ydim - averagingBlock, averagingBlock):
                cx = 0
                for x in range(0, xdim - averagingBlock, averagingBlock):
                    bE = ghostmap[y : y + averagingBlock, x : x + averagingBlock, c]
                    blkE[cy, cx, c] = np.mean(bE)
                    cx += 1
                cy += 1

        # normalize difference
        minval = np.min(blkE, axis=2)
        maxval = np.max(blkE, axis=2)
        for c in range(nQ):
            blkE[:, :, c] = (blkE[:, :, c] - minval) / (maxval - minval)

        # change plotsize (inches)
        self.plt.figure(figsize=(12, 8))

        if includeoriginal:
            sp = math.ceil(math.sqrt(nQ + 1))
            # Plot original image - needs to be normalized first for a subplot & converted to RGB, because cv2 works by default on BGR
            original_uint8 = cv2.convertScaleAbs(
                original
            )  # cv2.cvtColor expects input images to have depth of 8-bit per channel
            original_rgb = cv2.cvtColor(original_uint8, cv2.COLOR_BGR2RGB)
            originalRGB_normalized = original_rgb.astype(np.float32) / 255.0
            self.plt.subplot(sp, sp, 1)
            self.plt.imshow(originalRGB_normalized)
            self.plt.title("Original Image")
            self.plt.axis("off")
            # add maps:
            if grayscale:
                for c in range(nQ):
                    self.plt.subplot(sp, sp, c + 2)
                    self.plt.imshow(blkE[:, :, c], cmap="gray", vmin=0, vmax=1)
                    self.plt.axis("off")
                    self.plt.title("Quality " + str(Qmin + c * Qstep))
                    self.plt.draw()
            else:
                for c in range(nQ):
                    self.plt.subplot(sp, sp, c + 2)
                    self.plt.imshow(blkE[:, :, c], vmin=0, vmax=1)
                    self.plt.axis("off")
                    self.plt.title("Quality " + str(Qmin + c * Qstep))
                    self.plt.draw()
        else:
            sp = math.ceil(math.sqrt(nQ))
            if grayscale:
                for c in range(nQ):
                    self.plt.subplot(sp, sp, c + 1)
                    self.plt.imshow(blkE[:, :, c], cmap="gray", vmin=0, vmax=1)
                    self.plt.axis("off")
                    self.plt.title("Quality " + str(Qmin + c * Qstep))
                    self.plt.draw()
            else:
                for c in range(nQ):
                    self.plt.subplot(sp, sp, c + 1)
                    self.plt.imshow(blkE[:, :, c], vmin=0, vmax=1)
                    self.plt.axis("off")
                    self.plt.title("Quality " + str(Qmin + c * Qstep))
                    self.plt.draw()

        # Add main title
        self.plt.suptitle(
            "Ghost plots for grid offset X = "
            + str(shift_x)
            + " and Y = "
            + str(shift_y)
        )

        # Save plot directly to a file (temporarily)
        # Matplotlib deals with RGB images by default and does not provide functionality to switch between RGB and BGR
        # since sherloq primarily works with BGR, the conversion is made here using cv2 to read to image to a numpy array, BGR format, for further use
        temp_filename = "temp_ghostplot.png"
        plt.savefig(temp_filename, dpi=200)  # increase quality plot

        # Load the saved image in a numpy array, BGR format
        numpy_ghostplot = cv2.imread(temp_filename, cv2.IMREAD_COLOR)

        # save plot in memory so no recalculations are needed if user wants to revistit plot
        self.ghostmaps[shift_x + shift_y * 8] = [
            numpy_ghostplot,
            Qmin,
            Qmax,
            Qstep,
            includeoriginal,
            grayscale,
        ]

        # Remove the temporary file
        os.remove(temp_filename)

        # update viewer with plot:
        self.viewer.update_processed(numpy_ghostplot)

        self.plt.close()  # matplot figures are kept in memory unless closed
        self.process_button.setEnabled(True)  # allow new process to start
