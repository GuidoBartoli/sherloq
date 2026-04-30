from time import time
import numpy as np
import pickle
import cv2 as cv
from scipy.fftpack import dct
from scipy.ndimage import gaussian_filter
from skimage.transform import resize

from tools import ToolWidget
from viewer import ImageViewer
from PySide6.QtWidgets import QVBoxLayout, QLabel, QSlider, QHBoxLayout, QProgressBar
from PySide6.QtCore import Qt, QTimer

nameModel = r"DCT_model/RandomForest_2.0"
secSize = 16
rmFreq = 0

class DCTWidget(ToolWidget):
    def __init__(self, filename, image, parent=None):
        super(DCTWidget, self).__init__(parent)

        self.worker = None
        
        self.image = image

        with open(f"{nameModel}.pkl", "rb") as f:
            self.model = pickle.load(f)

        self.viewer = ImageViewer(self.image, self.image)
        
        self._pending_params = None

        self.label = QLabel("Initializing...")
        self.label_result = QLabel("")
        self.status_label = QLabel("")
        
        
        layout = QVBoxLayout()
        layout.addWidget(self.label)
        layout.addWidget(self.label_result)
        layout.addWidget(self.viewer)
        layout.addWidget(self.status_label)
        self.setLayout(layout)
        
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)

        layout.addWidget(self.progress_bar)
        
        self.update_timer = QTimer()
        self.update_timer.setInterval(200)
        self.update_timer.setSingleShot(True)
        self.update_timer.timeout.connect(self.process)
                
        # Slider Alpha
        self.alpha_slider = QSlider(Qt.Horizontal)
        self.alpha_slider.setRange(1, 100)
        self.alpha_slider.setValue(100) # default 1

        # Slide Threshold
        self.threshold_slider = QSlider(Qt.Horizontal)
        self.threshold_slider.setRange(10, 100) # 0.1 – 1
        self.threshold_slider.setValue(20) # default 0.2

        # Labels Sliders
        self.alpha_label = QLabel("Alpha: 1.0")
        self.threshold_label = QLabel("Threshold: 0.2")

        slider_layout = QHBoxLayout()
        slider_layout.addWidget(self.alpha_label)
        slider_layout.addWidget(self.alpha_slider)
        slider_layout.addWidget(self.threshold_label)
        slider_layout.addWidget(self.threshold_slider)

        layout.addLayout(slider_layout)
        
        
        self.live_timer = QTimer()
        self.live_timer.setInterval(100)
        self.live_timer.timeout.connect(self.update_live_time)
        
        self.time_label = QLabel("Time: -")
        layout.addWidget(self.time_label)
        
        self.alpha_slider.valueChanged.connect(self.update_parameters)
        self.threshold_slider.valueChanged.connect(self.update_parameters)

        self.preprocess()

    def update_parameters(self):
        self.alpha = self.alpha_slider.value() / 100
        self.threshold = self.threshold_slider.value() / 100

        self.alpha_label.setText(f"Alpha: {self.alpha:.2f}")
        self.threshold_label.setText(f"Threshold: {self.threshold:.2f}")

        self._pending_params = (self.alpha, self.threshold)
        self.update_timer.start()

    def update_progress(self, value):
        self.status_label.setText(f"Processing...")
        self.progress_bar.setValue(value)        
        
    def update_live_time(self):
        if hasattr(self, "start_time"):
            elapsed = time() - self.start_time
            progress = self.progress_bar.value()
            
            if progress < 5:
                self.time_label.setText(f"Time: {elapsed:.2f}s (estimating...)")
            
            elif progress > 0:
                eta = elapsed * (100 - progress) / progress

                if elapsed < 60:
                    self.time_label.setText(
                        f"Time: {elapsed:.2f}s (ETA: {eta:.2f}s)"
                    )
                else:
                    minutes = int(elapsed // 60)
                    seconds = elapsed % 60
                    self.time_label.setText(
                        f"Time: {minutes}m {seconds:.2f}s (ETA: {eta:.2f}s)"
                    )
            else:
                self.time_label.setText(f"Time: {elapsed:.2f}s")
            
    def on_finished(self, overlay, g, elapsed):
        self.viewer.update_processed(overlay)

        label = "Fake" if g >= 0.5 else "Real"
        self.label_result.setText(
            f"Classification: {label} ({g:.4f})"
        )

        self.progress_bar.setValue(100)
        self.status_label.setText("Done")

        minutes = int(elapsed // 60)
        seconds = elapsed % 60

        if elapsed < 60:
            self.time_label.setText(f"Time: {elapsed:.2f}s")
        else:
            minutes = int(elapsed // 60)
            seconds = elapsed % 60
            self.time_label.setText(f"Time: {minutes}m {seconds:.2f}s")
        
        if self._pending_params:
            self.alpha, self.threshold = self._pending_params
            self._pending_params = None
            
            QTimer.singleShot(0, self.process)
            
        self.live_timer.stop()
    
    def on_error(self, err):
        self.label.setText("Error occurred")
        self.progress_bar.setValue(0)

        print("DCT Error:", err)
    
    def extract_dct(self, img):
        img = dct(img, type=2, norm="ortho", axis=0)
        img = dct(img, type=2, norm="ortho", axis=1)

        img = np.abs(img)
        img += 1e-13
        img = np.log(img)

        img -= np.mean(img)
        img /= np.std(img)
        
        img = img[:secSize, :secSize]
        img[:rmFreq, :rmFreq] = 0

        feat_mean = np.mean(img)
        feat_std = np.std(img)
        feat_min = np.min(img)
        feat_max = np.max(img)

        return np.concatenate([
            img.flatten(),
            [feat_mean, feat_std, feat_min, feat_max]
        ])

    def predict(self):
        img_resized = resize(self.gray, (128, 128)).astype(np.float32)
        feat = self.extract_dct(img_resized)

        prob = self.model.predict_proba([feat])[0]
        g = prob[1]

        label = "Fake" if g >= 0.5 else "Real"

        self.label_result.setText(
            f"Classification: {label} ({g:.4f})"
        )

        return g

    def preprocess(self):
        self.label.setText("Preprocessing...")

        self.gray = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
        self.gray = self.gray.astype(np.float32) / 255.0
        self.gray = resize(self.gray, (256, 256)).astype(np.float32)
        QTimer.singleShot(0, self.process)

    def process(self):
        if self.worker and self.worker.isRunning():
            return
        
        self.label.setText("Model Output (Heatmap)")

        self.status_label.setText("Processing...")
        self.progress_bar.setValue(0)

        alpha = getattr(self, "alpha", 1.0)
        threshold = getattr(self, "threshold", 0.2)

        self.worker = DCTWorker(
            self.gray,
            self.image,
            self.model,
            self.extract_dct,
            alpha,
            threshold
        )

        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)

        self.start_time = time()
        self.time_label.setText("Time: 0.00s")

        self.live_timer.stop()
        self.live_timer.start()

        self.worker.start()
        
from PySide6.QtCore import QThread, Signal

class DCTWorker(QThread):
    progress = Signal(int)
    finished = Signal(object, float, float)
    error = Signal(str)

    def __init__(self, gray, image, model, extract_dct, alpha, threshold):
        super().__init__()
        self.gray = gray
        self.image = image
        self.model = model
        self.extract_dct = extract_dct
        self.alpha = alpha
        self.base_threshold = threshold

    def run(self):
        try:
            start = time()

            img_resized = resize(self.gray, (128, 128)).astype(np.float32)
            feat = self.extract_dct(img_resized)
            prob = self.model.predict_proba([feat])[0]
            g = prob[1]

            patch_size = 8
            stride = 8

            h, w = self.gray.shape
            heatmap = np.zeros((h, w))
            count = np.zeros((h, w))

            total = ((h - patch_size) // stride + 1) ** 2
            current = 0

            for y in range(0, h - patch_size + 1, stride):
                for x in range(0, w - patch_size + 1, stride):

                    patch = self.gray[y:y+patch_size, x:x+patch_size]
                    patch_resized = resize(patch, (128, 128)).astype(np.float32)

                    feat = self.extract_dct(patch_resized)
                    prob = self.model.predict_proba([feat])[0]

                    prob_fake = prob[1]
                    prob_fake *= (1 - self.alpha + self.alpha * g)

                    threshold = self.base_threshold - 0.2 * (g - 0.5)
                    prob_fake = max(0, prob_fake - threshold)
                    prob_fake = prob_fake ** 2

                    heatmap[y:y+patch_size, x:x+patch_size] += prob_fake
                    count[y:y+patch_size, x:x+patch_size] += 1

                    current += 1

                    if current % 10 == 0:
                        progress = int(current / total * 100)
                        self.progress.emit(progress)

            heatmap = (heatmap - heatmap.min()) / (heatmap.max() - heatmap.min() + 1e-8)
            heatmap = gaussian_filter(heatmap, sigma=3)
            heatmap = heatmap - heatmap.min()
            heatmap = heatmap / (heatmap.max() + 1e-8)

            heatmap = (heatmap * 255).astype(np.uint8)
            heatmap_color = cv.applyColorMap(heatmap, cv.COLORMAP_JET)

            heatmap_color = cv.resize(
                heatmap_color,
                (self.image.shape[1], self.image.shape[0])
            )

            overlay = cv.addWeighted(
                self.image.astype(np.uint8),
                0.6,
                heatmap_color.astype(np.uint8),
                0.4,
                0
            )

            elapsed = time() - start

            self.finished.emit(overlay, g, elapsed)

        except Exception as e:
            self.error.emit(str(e))