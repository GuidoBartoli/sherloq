# This code implements probability maps and fourier maps for detecting traces of resampling as explained in the paper:
# "Exposing Digital Forgeries by Detecting Traces of Re-sampling" by Hany Farid & Alin C. Popescu
# The book "Photo Forensics" by Hany Farid gives a more detailed explanation of the technique for those interested

from PySide6.QtWidgets import (
    QWidget,
    QDoubleSpinBox,
    QVBoxLayout,
    QSlider,
    QHBoxLayout,
    QLabel,
    QSpinBox,
    QCheckBox,
    QPushButton,
    QGridLayout,
    QSizePolicy,
    QScrollArea,
)
from PySide6.QtGui import QTransform
from PySide6.QtCore import Qt
from tools import ToolWidget
from viewer import ImageViewer

# resampling nessecary imports
import numpy as np
import cv2
import random
import matplotlib.pyplot as plt

plt.rcParams["savefig.dpi"] = 200  # increase deafault save resolution
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.backend_bases import MouseButton
from matplotlib.figure import Figure
from utility import modify_font


class ResamplingWidget(ToolWidget):
    # tool layout
    def __init__(self, filename, image, parent=None):
        super(ResamplingWidget, self).__init__(parent)

        self.filename = filename
        self.image_original = image
        # read and normalize image into grayscale for analysis
        self.imagegray = cv2.imread(filename, cv2.IMREAD_GRAYSCALE)
        if self.imagegray is None:
            error_label = QLabel(self.tr("Unable to detect resampling on raw images!"))
            modify_font(error_label, bold=True)
            error_label.setStyleSheet("color: #FF0000")
            error_label.setAlignment(Qt.AlignCenter)
            main_layout = QVBoxLayout()
            main_layout.addWidget(error_label)
            self.setLayout(main_layout)
            return
        self.imagegray = self.imagegray - self.imagegray.min()
        self.imagegray = self.imagegray / self.imagegray.max()

        self.imagegray_nomalized_copy = self.imagegray.copy()
        self.imagegray_copy_for_probabilitymaps = self.imagegray.copy()
        self.userfourplot = self.imagegray.copy()

        self.selected_points_probability = []
        self.selected_points_fourier = []
        self.probability_maps = []
        self.fourier_maps = []
        self.original_sizes_widgets = {}

        # prepare user interface
        self.calculate_probability_button = QPushButton(
            self.tr("Calculate probability")
        )
        self.calculate_fourier_button = QPushButton(self.tr("Calculate fourier"))

        # N x N symmetric interpolation model to estimate interpolation probability
        self.filter_3x3_Check = QCheckBox(self.tr("3x3 probability filter"))
        self.filter_3x3_Check.setChecked(True)
        self.filter_5x5_Check = QCheckBox(self.tr("5x5 probability filter"))
        self.filter_5x5_Check.setChecked(False)

        # bind toggle
        self.filter_3x3_Check.clicked.connect(
            lambda: self.filter_5x5_Check.setChecked(False)
        )
        self.filter_5x5_Check.clicked.connect(
            lambda: self.filter_3x3_Check.setChecked(False)
        )

        # pointpicker probability
        self.probability_check = QCheckBox(self.tr("Probability Windows"))
        self.probability_check.setChecked(True)

        # pointpicker fourier
        self.fourier_check = QCheckBox(self.tr("Fourier Windows"))
        self.fourier_check.setChecked(False)

        # combine top layout
        top_layout = QGridLayout()
        top_layout.addWidget(
            QLabel(
                self.tr(
                    "A resampling analysis is computationally expensive and for a\n"
                    + "high resoluion image (1280x720 pixels), it can take a few minutes to process.\n"
                    + "It is more efficient and more effective to calculate the probability map for only the suspected area.\n"
                    + "This is because a small resampled area is more readily "
                    + "detected when the larger image is not considered by the algoritm."
                )
            ),
            0,
            0,
        )

        top_layout.addWidget(
            QLabel(
                self.tr(
                    "A 5x5 probability filter might be able to uncover more complex interpolation algorithms, but the computing time for the probability maps increases significantly.\n"
                    + "Left mouse button: choose location to construct area of interest\n"
                    + "Right mouse button: delete last location."
                )
            ),
            1,
            0,
        )

        top_layout.addWidget(
            QLabel(
                self.tr(
                    "If no area of interest is chosen for probability, the probability map for the entire image will be calculated.\n"
                    + "Fourier maps will automatically be calculated for each area of interest by the probability map + for smaller sub windows applied when 'Fourier Windows' is cheched.\n"
                    "Please make sure areas do not overlap for probability maps!"
                )
            ),
            2,
            0,
        )

        top_layout.addWidget(self.probability_check, 0, 1)
        top_layout.addWidget(self.fourier_check, 0, 2)

        checkbox_layout = QVBoxLayout()
        checkbox_layout.addWidget(self.filter_3x3_Check)
        checkbox_layout.addWidget(self.filter_5x5_Check)
        top_layout.addLayout(checkbox_layout, 1, 1)

        top_layout.addWidget(self.calculate_probability_button, 2, 1)
        top_layout.addWidget(self.calculate_fourier_button, 2, 2)

        # fourier maps contruction parameters
        four_parameters_layout = QGridLayout()
        # pre-processing option: hanning or rotationally invariant window?
        self.hanning_check = QCheckBox(self.tr("Hanning"))
        self.hanning_check.setChecked(True)
        self.rotationally_invariant_window_check = QCheckBox(self.tr("R. I. W."))
        self.rotationally_invariant_window_check.setChecked(False)
        # bind toggle options
        self.hanning_check.clicked.connect(
            lambda: self.rotationally_invariant_window_check.setChecked(False)
        )
        self.rotationally_invariant_window_check.clicked.connect(
            lambda: self.hanning_check.setChecked(False)
        )

        # upsample and take center options
        self.upsample_check = QCheckBox(self.tr("Upsample"))
        self.upsample_check.setChecked(True)
        self.center_four_check = QCheckBox(self.tr("Take center of Fourier"))
        self.center_four_check.setChecked(False)

        # filter options
        self.simple_highpass_check = QCheckBox(self.tr("Highpass 1"))
        self.simple_highpass_check.setChecked(True)
        self.complex_highpass_check = QCheckBox(self.tr("Highpass 2"))
        self.complex_highpass_check.setChecked(False)
        self.simple_highpass_check.clicked.connect(
            lambda: self.complex_highpass_check.setChecked(False)
        )
        self.complex_highpass_check.clicked.connect(
            lambda: self.simple_highpass_check.setChecked(False)
        )

        # gamma correction
        self.gamma_spin = QDoubleSpinBox()
        self.gamma_spin.setRange(0, 5)
        self.gamma_spin.setSingleStep(0.1)
        self.gamma_spin.setValue(4)

        # rescale option
        self.rescale_check = QCheckBox(self.tr("Rescale Spectrum"))
        self.rescale_check.setChecked(True)

        # assemble four parameters
        four_parameters_layout.addWidget(self.hanning_check, 0, 0)
        four_parameters_layout.addWidget(self.rotationally_invariant_window_check, 1, 0)
        four_parameters_layout.addWidget(self.upsample_check, 1, 1)
        four_parameters_layout.addWidget(self.center_four_check, 1, 2)
        four_parameters_layout.addWidget(self.simple_highpass_check, 0, 3)
        four_parameters_layout.addWidget(self.complex_highpass_check, 1, 3)
        four_parameters_layout.addWidget(QLabel(self.tr("Gamma correction")), 0, 4)
        four_parameters_layout.addWidget(self.gamma_spin, 1, 4)
        four_parameters_layout.addWidget(self.rescale_check, 1, 5)

        # canvas with picture for probability
        figure_prob = Figure()
        self.canvas = FigureCanvas(figure_prob)
        self.axes = self.canvas.figure.subplots()
        # imshow creates an object, so the object is saved here so it can be updated in the future to show new content
        # repeatedly using imshow without clearing data will result in increased memory and subsequently slower performance
        self.probability_image_canvas_object = self.axes.imshow(
            self.imagegray, cmap="gray"
        )
        self.canvas.mpl_connect("button_press_event", self.click_on_canvas)
        self.canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.axes.figure.canvas.draw()

        # add toolbar to probability figure for soom and export
        self.toolbar_prob = NavigationToolbar(self.canvas, self)

        # canvas for fourier_Maps
        self.figure_four = plt.figure(figsize=(10, 15))
        self.canvas_fourier_maps = FigureCanvas(self.figure_four)
        self.canvas_fourier_maps.mpl_connect("button_press_event", self.click_on_canvas)
        self.canvas_fourier_maps.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        # add toolbar to fourier figure for zoom and export
        self.toolbar_four = NavigationToolbar(self.canvas_fourier_maps, self)

        # bind user interactions
        self.calculate_probability_button.clicked.connect(
            self.calculate_probability_maps
        )
        self.calculate_fourier_button.clicked.connect(self.calculate_fourier_maps)
        self.probability_check.clicked.connect(
            lambda: self.fourier_check.setChecked(False)
        )
        self.fourier_check.clicked.connect(
            lambda: self.probability_check.setChecked(False)
        )

        # assemble layout
        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        self.four_parameters_label = QLabel(
            "Fourier map options: (Make sure to have Hanning or R. I. W. checked, same for Highpass 1 or 2.)"
        )
        main_layout.addWidget(self.four_parameters_label)
        main_layout.addLayout(four_parameters_layout)

        # add slider for zooming
        self.zoom_label = QLabel("Zoom:")
        self.zoom_slider = QSlider(Qt.Horizontal)
        self.zoom_slider.setMinimum(100)
        self.zoom_slider.setMaximum(200)
        self.zoom_slider.setValue(100)
        self.zoom_slider.setTickInterval(5)
        self.zoom_slider.valueChanged.connect(self.slider_zoom)

        main_layout.addWidget(self.zoom_label)
        main_layout.addWidget(self.zoom_slider)

        # make the figures scrollable
        self.scroll_area = QScrollArea(self)
        self.scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(self.scroll_widget)
        self.scroll_area.setWidget(self.scroll_widget)
        self.scroll_area.setWidgetResizable(True)

        scroll_layout.addWidget(self.toolbar_prob)
        scroll_layout.addWidget(self.canvas)
        self.original_sizes_widgets[self.canvas] = self.canvas.sizeHint()

        scroll_layout.addWidget(self.toolbar_four)
        scroll_layout.addWidget(self.canvas_fourier_maps)
        self.original_sizes_widgets[
            self.canvas_fourier_maps
        ] = self.canvas_fourier_maps.sizeHint()

        main_layout.addWidget(self.scroll_area)

        self.setLayout(main_layout)

    def slider_zoom(self):
        zoom_factor = self.zoom_slider.value() / 100
        for child_widget in self.scroll_widget.findChildren(QWidget):
            if child_widget in self.original_sizes_widgets:
                original_size = self.original_sizes_widgets[child_widget]
                new_width = original_size.width() * zoom_factor
                new_height = original_size.height() * zoom_factor
                child_widget.setFixedSize(new_width, new_height)

    def calculate_probability_maps(self):
        # to do, if probability maps already calculated, do not calculate again;
        # untill then, make fourier maps empty before calculating again:
        self.probability_maps = []
        self.imagegray_copy_for_probabilitymaps = self.imagegray_nomalized_copy.copy()
        self.calculate_probability_button.setEnabled(False)  # wait for processing

        if len(self.selected_points_probability) == 0:
            if self.filter_5x5_Check.isChecked():
                processed_part = self.calculate_probability_map_5x5(
                    self.imagegray_nomalized_copy
                )
            else:
                processed_part = self.calculate_probability_map_3x3(
                    self.imagegray_nomalized_copy
                )

            self.probability_maps.append(processed_part)
            self.imagegray_copy_for_probabilitymaps = processed_part
        else:
            for i in range(0, len(self.selected_points_probability), 2):
                if i + 1 < len(self.selected_points_probability):
                    (x1, y1) = self.selected_points_probability[i]
                    (x2, y2) = self.selected_points_probability[i + 1]
                    if self.filter_5x5_Check.isChecked():
                        processed_part = self.calculate_probability_map_5x5(
                            self.imagegray_nomalized_copy[
                                min(y1, y2) : max(y1, y2) + 1,
                                min(x1, x2) : max(x1, x2) + 1,
                            ]
                        )
                        self.probability_maps.append(processed_part)
                        self.imagegray_copy_for_probabilitymaps[
                            min(y1, y2) + 2 : max(y1, y2) - 1,
                            min(x1, x2) + 2 : max(x1, x2) - 1,
                        ] = processed_part
                    else:
                        processed_part = self.calculate_probability_map_3x3(
                            self.imagegray_nomalized_copy[
                                min(y1, y2) : max(y1, y2) + 1,
                                min(x1, x2) : max(x1, x2) + 1,
                            ]
                        )
                        self.probability_maps.append(processed_part)
                        self.imagegray_copy_for_probabilitymaps[
                            min(y1, y2) + 1 : max(y1, y2), min(x1, x2) + 1 : max(x1, x2)
                        ] = processed_part

        # update content instead of using imshow again
        self.probability_image_canvas_object.set_data(
            self.imagegray_copy_for_probabilitymaps
        )
        self.axes.figure.canvas.draw()
        self.calculate_probability_button.setEnabled(True)

    def calculate_probability_map_3x3(self, process_part):
        a = np.random.rand(8)
        a = a / a.sum()
        s = 0.005
        d = 0.1
        F, f = self.build_matrices_for_processing_3x3(process_part)
        c = 0
        w = np.zeros(F.shape[0])

        while c < 100:
            s2 = 0
            k = 0

            for y in range(1, process_part.shape[0] - 1):
                for x in range(1, process_part.shape[1] - 1):
                    r = self.compute_residual_3x3(a, process_part, x, y)
                    g = np.exp(-(r ** 2) / s)
                    w[k] = g / (g + d)
                    s2 = s2 + w[k] * r ** 2
                    k = k + 1

            s = s2 / w.sum()
            a2 = np.linalg.inv(F.T * w * w @ F) @ F.T * w * w @ f

            if np.linalg.norm(a - a2) < 0.01:
                break
            else:
                a = a2
                c = c + 1

        return w.reshape(process_part.shape[0] - 2, process_part.shape[1] - 2)

    def calculate_probability_map_5x5(self, process_part):
        a = np.random.rand(24)
        a = a / a.sum()
        s = 0.005
        d = 0.1
        F, f = self.build_matrices_for_processing_5x5(process_part)
        c = 0
        w = np.zeros(F.shape[0])

        while c < 100:
            s2 = 0
            k = 0

            for y in range(2, process_part.shape[0] - 2):
                for x in range(2, process_part.shape[1] - 2):
                    r = self.compute_residual_5x5(a, process_part, x, y)
                    g = np.exp(-(r ** 2) / s)
                    w[k] = g / (g + d)
                    s2 = s2 + w[k] * r ** 2
                    k = k + 1

            s = s2 / w.sum()
            a2 = np.linalg.inv(F.T * w * w @ F) @ F.T * w * w @ f

            if np.linalg.norm(a - a2) < 0.01:
                break
            else:
                a = a2
                c = c + 1

        return w.reshape(process_part.shape[0] - 4, process_part.shape[1] - 4)

    def calculate_fourier_maps(self):
        # to do, if fouriermap already calculated, do not calculate again
        # untill then, make fourier maps empty before calculating again:
        self.fourier_maps = []
        self.calculate_fourier_button.setEnabled(False)  # wait for processing

        num_plots = int(len(self.selected_points_fourier) / 2) + len(
            self.probability_maps
        )
        self.canvas_fourier_maps.figure.clf()
        if num_plots > 1:
            self.axes_fourier_maps = self.canvas_fourier_maps.figure.subplots(
                num_plots, 2
            )
        else:
            self.axes_fourier_maps = np.array(
                [
                    [
                        self.canvas_fourier_maps.figure.add_subplot(1, 2, 1),
                        self.canvas_fourier_maps.figure.add_subplot(1, 2, 2),
                    ]
                ]
            )

        i = 0
        counter_selected_prob_points = 0
        for prob_map in self.probability_maps:
            if (
                len(self.probability_maps) == 1
                and len(self.selected_points_probability) == 0
            ):
                self.fourier_maps.append(self.calculate_fourier_map(prob_map))
                self.axes_fourier_maps[i, 0].imshow(
                    prob_map, cmap="gray", vmin=0, vmax=1
                )
                self.axes_fourier_maps[i, 0].axis("off")
                self.axes_fourier_maps[i, 1].imshow(
                    self.fourier_maps[i], cmap="gray", vmin=0, vmax=1
                )
                self.axes_fourier_maps[i, 1].axis("off")
                i = i + 1
            else:
                (x1, y1) = self.selected_points_probability[
                    counter_selected_prob_points
                ]
                (x2, y2) = self.selected_points_probability[
                    counter_selected_prob_points + 1
                ]
                self.fourier_maps.append(self.calculate_fourier_map(prob_map))
                self.userfourplot = self.imagegray_copy_for_probabilitymaps.copy()
                self.userfourplot[
                    min(y1, y2) : max(y1, y2) + 1,
                    [
                        x1,
                        (x1 + 1),
                        (x1 - 1),
                        (x1 - 2),
                        (x1 + 2),
                        x2,
                        (x2 + 1),
                        (x2 - 1),
                        (x2 - 2),
                        (x2 + 2),
                    ],
                ] = 0  # y-stripes
                self.userfourplot[
                    [
                        y1,
                        (y1 - 1),
                        (y1 + 1),
                        (y1 - 2),
                        (y1 + 2),
                        y2,
                        (y2 - 1),
                        (y2 + 1),
                        (y2 - 2),
                        (y2 + 2),
                    ],
                    min(x1, x2) : max(x1, x2) + 1,
                ] = 0  # x-stripes
                self.axes_fourier_maps[i, 0].imshow(
                    self.userfourplot, cmap="gray", vmin=0, vmax=1
                )
                self.axes_fourier_maps[i, 0].axis("off")
                self.axes_fourier_maps[i, 1].imshow(
                    self.fourier_maps[i], cmap="gray", vmin=0, vmax=1
                )
                self.axes_fourier_maps[i, 1].axis("off")
                i = i + 1
                counter_selected_prob_points = counter_selected_prob_points + 2
        for j in range(0, len(self.selected_points_fourier), 2):
            if j + 1 < len(self.selected_points_fourier):
                (x1, y1) = self.selected_points_fourier[j]
                (x2, y2) = self.selected_points_fourier[j + 1]

                processed_part = self.calculate_fourier_map(
                    self.imagegray_copy_for_probabilitymaps[
                        min(y1, y2) : max(y1, y2) + 1, min(x1, x2) : max(x1, x2) + 1
                    ]
                )
                self.fourier_maps.append(processed_part)
                self.userfourplot = self.imagegray_copy_for_probabilitymaps.copy()
                self.userfourplot[
                    min(y1, y2) : max(y1, y2) + 1,
                    [
                        x1,
                        (x1 + 1),
                        (x1 - 1),
                        (x1 - 2),
                        (x1 + 2),
                        x2,
                        (x2 + 1),
                        (x2 - 1),
                        (x2 - 2),
                        (x2 + 2),
                    ],
                ] = 0  # y-stripes
                self.userfourplot[
                    [
                        y1,
                        (y1 - 1),
                        (y1 + 1),
                        (y1 - 2),
                        (y1 + 2),
                        y2,
                        (y2 - 1),
                        (y2 + 1),
                        (y2 - 2),
                        (y2 + 2),
                    ],
                    min(x1, x2) : max(x1, x2) + 1,
                ] = 0  # x-stripes
                self.axes_fourier_maps[i, 0].imshow(
                    self.userfourplot, cmap="gray", vmin=0, vmax=1
                )
                self.axes_fourier_maps[i, 0].axis("off")
                self.axes_fourier_maps[i, 1].imshow(
                    self.fourier_maps[i], cmap="gray", vmin=0, vmax=1
                )
                self.axes_fourier_maps[i, 1].axis("off")
                i = i + 1

        self.figure_four.subplots_adjust(
            left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0.05, hspace=0.05
        )

        self.canvas_fourier_maps.figure.canvas.draw()
        self.calculate_fourier_button.setEnabled(True)

    def make_rotational_invariant_window(self, shape):
        rows, cols = shape
        W = np.zeros((rows, cols))
        center_x, center_y = rows // 2, cols // 2
        max_radius = np.sqrt(center_x ** 2 + center_y ** 2)
        for i in range(rows):
            for j in range(cols):
                r = (
                    np.sqrt((i - center_x) ** 2 + (j - center_y) ** 2)
                    / max_radius
                    * np.sqrt(2)
                )
                if r < 3 / 4:
                    W[i, j] = 1
                elif r <= np.sqrt(2):
                    W[i, j] = 0.5 + 0.5 * np.cos(
                        np.pi * (r - 3 / 4) / (np.sqrt(2) - 3 / 4)
                    )

        return W

    def make_high_pass_filter(self, shape):
        rows, cols = shape
        H = np.zeros((rows, cols))
        center_x, center_y = rows // 2, cols // 2
        max_radius = np.sqrt(center_x ** 2 + center_y ** 2)
        for i in range(rows):
            for j in range(cols):
                r = (
                    np.sqrt((i - center_x) ** 2 + (j - center_y) ** 2)
                    / max_radius
                    * np.sqrt(2)
                )
                if r <= np.sqrt(2):
                    H[i, j] = 0.5 - 0.5 * np.cos(np.pi * r / np.sqrt(2))

        return H

    def calculate_fourier_map(self, process_prob_map):
        # take center square portion of probability map, this will make analysis easier and more consistent
        x, y = process_prob_map.shape
        size = min(x, y)
        half_size = size // 2
        center_x, center_y = x // 2, y // 2
        square_prob_map = process_prob_map[
            center_x - half_size : center_x + half_size,
            center_y - half_size : center_y + half_size,
        ]

        # apply pre-processing window
        if self.hanning_check.isChecked():
            # Apply Hanning window
            hanning_window = np.hanning(square_prob_map.shape[0])[:, None] * np.hanning(
                square_prob_map.shape[1]
            )
            windowed_prob_map = square_prob_map * hanning_window
        elif self.rotationally_invariant_window_check.isChecked():
            # apply rot. invariant window
            W = self.make_rotational_invariant_window(square_prob_map.shape)
            windowed_prob_map = square_prob_map * W
        else:
            # if none is checked, return
            return

        if self.upsample_check.isChecked():
            upsampled = cv2.pyrUp(windowed_prob_map)
        else:
            upsampled = windowed_prob_map

        # fourier transform
        dft = np.fft.fft2(upsampled)
        fourier = np.fft.fftshift(dft)

        # take center?
        if self.center_four_check.isChecked():
            height, width = fourier.shape
            center_x, center_y = width // 2, height // 2
            half_size = width // 4
            fourier_center = fourier[
                center_y - half_size : center_y + half_size,
                center_x - half_size : center_x + half_size,
            ]
        else:
            fourier_center = fourier

        # filter option
        if self.simple_highpass_check.isChecked():
            # create circular mask:
            rows, cols = fourier_center.shape
            center = (int(cols / 2), int(rows / 2))
            radius = int(0.1 * (min(rows, cols) / 2))
            if radius == 0:
                radius = 1  # radius should at least be = 1

            Y, X = np.ogrid[:rows, :cols]
            dist_from_center = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2)

            mask = dist_from_center <= radius
            filtered_spectrum = fourier_center.copy()
            filtered_spectrum[mask] = 0

        elif self.complex_highpass_check.isChecked():
            H = self.make_high_pass_filter(fourier_center.shape)
            filtered_spectrum = fourier_center * H
        else:
            return

        # scale and gamma correct
        magnitude_spectrum = np.abs(filtered_spectrum)
        scaled_spectrum = (magnitude_spectrum - magnitude_spectrum.min()) / (
            magnitude_spectrum.max() - magnitude_spectrum.min()
        )

        gamma_corrected_spectrum = np.power(scaled_spectrum, self.gamma_spin.value())

        # rescale if checked
        if self.rescale_check.isChecked():
            rescaled_spectrum = gamma_corrected_spectrum * magnitude_spectrum.max()
        else:
            rescaled_spectrum = gamma_corrected_spectrum

        return rescaled_spectrum

    def click_on_canvas(self, event):
        if self.toolbar_four.mode == "" and self.toolbar_prob.mode == "":
            if event.inaxes:
                if self.probability_check.isChecked():
                    if event.button == MouseButton.LEFT:
                        x, y = int(event.xdata), int(event.ydata)
                        self.selected_points_probability.append((x, y))
                        print(f"Selected point: ({x}, {y})")
                        self.imagegray = self.imagegray_nomalized_copy.copy()
                        self.draw_selection(self.selected_points_probability)
                        self.probability_image_canvas_object.set_data(self.imagegray)
                        self.axes.figure.canvas.draw()
                    elif (
                        event.button == MouseButton.RIGHT
                        and self.selected_points_probability
                    ):
                        print(self.selected_points_probability)
                        self.selected_points_probability.pop()
                        print(self.selected_points_probability)
                        self.imagegray = self.imagegray_nomalized_copy.copy()
                        self.draw_selection(self.selected_points_probability)
                        self.probability_image_canvas_object.set_data(self.imagegray)
                        self.axes.figure.canvas.draw()
                elif self.fourier_check.isChecked():
                    if event.button == MouseButton.LEFT:
                        x, y = int(event.xdata), int(event.ydata)
                        self.selected_points_fourier.append((x, y))
                        self.userfourplot = (
                            self.imagegray_copy_for_probabilitymaps.copy()
                        )
                        self.draw_selection(self.selected_points_fourier)
                        self.probability_image_canvas_object.set_data(self.userfourplot)
                        self.axes.figure.canvas.draw()
                    elif (
                        event.button == MouseButton.RIGHT
                        and self.selected_points_fourier
                    ):
                        self.selected_points_fourier.pop()
                        self.userfourplot = (
                            self.imagegray_copy_for_probabilitymaps.copy()
                        )
                        self.draw_selection(self.selected_points_fourier)
                        self.probability_image_canvas_object.set_data(self.userfourplot)
                        self.axes.figure.canvas.draw()

    def draw_selection(self, selected_points):
        i = 1
        if self.probability_check.isChecked():
            for point in selected_points:
                cv2.circle(self.imagegray, point, 3, (1,), -1)
                if (i % 2) == 0:
                    self.draw_rectangle(selected_points[i - 2], selected_points[i - 1])

                i = i + 1
        elif self.fourier_check.isChecked():
            for point in selected_points:
                cv2.circle(self.userfourplot, point, 3, (1,), -1)
                if (i % 2) == 0:
                    self.draw_rectangle(selected_points[i - 2], selected_points[i - 1])

                i = i + 1

    def draw_rectangle(self, p1, p2):
        (x1, y1) = p1
        (x2, y2) = p2
        if self.probability_check.isChecked():
            self.imagegray[
                min(y1, y2) : max(y1, y2) + 1,
                [x1, (x1 + 1), (x1 - 1), x2, (x2 + 1), (x2 - 1)],
            ] = 0  # y-stripes
            self.imagegray[
                [y1, (y1 - 1), (y1 + 1), y2, (y2 - 1), (y2 + 1)],
                min(x1, x2) : max(x1, x2) + 1,
            ] = 0  # x-stripes
        elif self.fourier_check.isChecked():
            self.userfourplot[
                min(y1, y2) : max(y1, y2) + 1,
                [x1, (x1 + 1), (x1 - 1), x2, (x2 + 1), (x2 - 1)],
            ] = 0  # y-stripes
            self.userfourplot[
                [y1, (y1 - 1), (y1 + 1), y2, (y2 - 1), (y2 + 1)],
                min(x1, x2) : max(x1, x2) + 1,
            ] = 0  # x-stripes

    def build_matrices_for_processing_3x3(self, I):
        k = 0
        F = np.zeros(
            ((I.shape[0] - 2) * (I.shape[1] - 2), 8)
        )  # the i^th row of the matrix F corresponds to the pixel neighbors of the i^th element of f
        f = np.zeros(
            (I.shape[0] - 2) * (I.shape[1] - 2)
        )  # the vector f is the image strung out in row-order
        for y in range(1, I.shape[0] - 1):
            for x in range(1, I.shape[1] - 1):
                F[k, 0] = I[y - 1, x - 1]
                F[k, 1] = I[y - 1, x]
                F[k, 2] = I[y - 1, x + 1]
                F[k, 3] = I[y, x - 1]
                # skip I(x; y) corresponding to a(3; 3) which is 0
                F[k, 4] = I[y, x + 1]
                F[k, 5] = I[y + 1, x - 1]
                F[k, 6] = I[y + 1, x]
                F[k, 7] = I[y + 1, x + 1]
                f[k] = I[y, x]
                k = k + 1
        return F, f

    def build_matrices_for_processing_5x5(self, I):
        k = 0
        F = np.zeros(
            ((I.shape[0] - 4) * (I.shape[1] - 4), 24)
        )  # the i^th row of the matrix F corresponds to the pixel neighbors of the i^th element of f
        f = np.zeros(
            (I.shape[0] - 4) * (I.shape[1] - 4)
        )  # the vector f is the image strung out in row-order
        for y in range(2, I.shape[0] - 2):
            for x in range(2, I.shape[1] - 2):
                F[k, 0] = I[y - 2, x - 2]
                F[k, 1] = I[y - 2, x - 1]
                F[k, 2] = I[y - 2, x]
                F[k, 3] = I[y - 2, x + 1]
                F[k, 4] = I[y - 2, x + 2]
                F[k, 5] = I[y - 1, x - 2]
                F[k, 6] = I[y - 1, x - 1]
                F[k, 7] = I[y - 1, x]
                F[k, 8] = I[y - 1, x + 1]
                F[k, 9] = I[y - 1, x + 2]
                F[k, 10] = I[y, x - 2]
                F[k, 11] = I[y, x - 1]
                # skip I(x; y) corresponding to a(3; 3) which is 0
                F[k, 12] = I[y, x + 1]
                F[k, 13] = I[y, x + 2]
                F[k, 14] = I[y + 1, x - 2]
                F[k, 15] = I[y + 1, x - 1]
                F[k, 16] = I[y + 1, x]
                F[k, 17] = I[y + 1, x + 1]
                F[k, 18] = I[y + 1, x + 2]
                F[k, 19] = I[y + 2, x - 2]
                F[k, 20] = I[y + 2, x - 1]
                F[k, 21] = I[y + 2, x]
                F[k, 22] = I[y + 2, x + 1]
                F[k, 23] = I[y + 2, x + 2]
                f[k] = I[y, x]
                k = k + 1
        return F, f

    def compute_residual_3x3(self, a, I, x, y):
        r = I[y, x] - (
            a[0] * I[y - 1, x - 1]
            + a[1] * I[y - 1, x]
            + a[2] * I[y - 1, x + 1]
            + a[3] * I[y, x - 1]
            + a[4] * I[y, x + 1]
            + a[5] * I[y + 1, x - 1]
            + a[6] * I[y + 1, x]
            + a[7] * I[y + 1, x + 1]
        )
        return r

    def compute_residual_5x5(self, a, I, x, y):
        r = I[y, x] - (
            a[0] * I[y - 2, x - 2]
            + a[1] * I[y - 2, x - 1]
            + a[2] * I[y - 2, x]
            + a[3] * I[y - 2, x + 1]
            + a[4] * I[y - 2, x + 2]
            + a[5] * I[y - 1, x - 2]
            + a[6] * I[y - 1, x - 1]
            + a[7] * I[y - 1, x]
            + a[8] * I[y - 1, x + 1]
            + a[9] * I[y - 1, x + 2]
            + a[10] * I[y, x - 2]
            + a[11] * I[y, x - 1]
            + a[12] * I[y, x + 1]
            + a[13] * I[y, x + 2]
            + a[14] * I[y + 1, x - 2]
            + a[15] * I[y + 1, x - 1]
            + a[16] * I[y + 1, x]
            + a[17] * I[y + 1, x + 1]
            + a[18] * I[y + 1, x + 2]
            + a[19] * I[y + 2, x - 2]
            + a[20] * I[y + 2, x - 1]
            + a[21] * I[y + 2, x]
            + a[22] * I[y + 2, x + 1]
            + a[23] * I[y + 2, x + 2]
        )
        return r
