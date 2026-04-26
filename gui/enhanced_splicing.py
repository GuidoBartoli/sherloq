import cv2 as cv
import numpy as np
from time import time
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (QGridLayout, QHBoxLayout, QVBoxLayout, QPushButton, 
                               QMessageBox, QTabWidget, QSlider, QLabel, QWidget)

from skimage.segmentation import slic, mark_boundaries
from tools import ToolWidget
from viewer import ImageViewer
from utility import modify_font

# --- 1. THE MATH WORKER THREAD ---
class ConsensusWorker(QThread):
    finished = Signal(object, object, object, object, object, object)
    error = Signal(str)

    def __init__(self, image_array):
        super().__init__()
        self.image_array = image_array

    def run(self):
        try:
            img_bgr = self.image_array
            img_rgb = cv.cvtColor(img_bgr, cv.COLOR_BGR2RGB)
            img_gray = cv.cvtColor(img_bgr, cv.COLOR_BGR2GRAY).astype(np.float64)
            
            # --- NEW: Convert to YCbCr to isolate pure color! ---
            img_ycbcr = cv.cvtColor(img_bgr, cv.COLOR_BGR2YCrCb)
            # Split out the Cr (Red-difference) and Cb (Blue-difference) channels
            _, cr_channel, cb_channel = cv.split(img_ycbcr)
            
            # Create a combined color variance map
            color_variance_map = np.abs(cr_channel.astype(np.float64) - cb_channel.astype(np.float64))
            
            # 1. SLIC Segmentation
            segments = slic(img_rgb, n_segments=999, compactness=10, sigma=1, start_label=1)
            segment_ids = np.unique(segments)

            # 2. Pre-computing Global Math
            noise_matrix = cv.Laplacian(img_gray, cv.CV_64F)
            
            diff_h = np.abs(img_gray[:, :-1] - img_gray[:, 1:])
            diff_v = np.abs(img_gray[:-1, :] - img_gray[1:, :])
            grid_map = np.zeros_like(img_gray)
            for i in range(7, img_gray.shape[0] - 1, 8):
                grid_map[i, :] += diff_v[i, :]
            for j in range(7, img_gray.shape[1] - 1, 8):
                grid_map[:, j] += diff_h[:, j]
                
            # Compress the image in RAM at 90% quality
            _, encoded_img = cv.imencode('.jpg', img_bgr, [int(cv.IMWRITE_JPEG_QUALITY), 90])
            compressed_img = cv.imdecode(encoded_img, 1)
            
            # Calculate the absolute difference and convert to grayscale
            ela_diff = cv.absdiff(img_bgr, compressed_img)
            ela_gray = cv.cvtColor(ela_diff, cv.COLOR_BGR2GRAY).astype(np.float64)
            
            # 3. Extracting Features
            noise_heatmap = np.zeros_like(img_gray)
            dct_heatmap = np.zeros_like(img_gray)
            freq_heatmap = np.zeros_like(img_gray)
            ela_heatmap = np.zeros_like(img_gray)
            color_heatmap = np.zeros_like(img_gray)
            
            for seg_id in segment_ids:
                mask = (segments == seg_id)
                noise_heatmap[mask] = np.var(noise_matrix[mask])
                dct_heatmap[mask] = np.mean(grid_map[mask])
                ela_heatmap[mask] = np.mean(ela_gray[mask])
                color_heatmap[mask] = np.var(color_variance_map[mask])
                 
                y, x = np.where(mask)
                if len(y) > 0 and len(x) > 0:
                    roi = img_gray[np.min(y):np.max(y)+1, np.min(x):np.max(x)+1]
                    if roi.shape[0] > 1 and roi.shape[1] > 1:
                        f_transform = np.fft.fft2(roi)
                        f_shift = np.fft.fftshift(f_transform)
                        freq_heatmap[mask] = np.var(np.abs(f_shift))

            def normalize_heatmap(hm):
                hm_min, hm_max = hm.min(), hm.max()
                if hm_max - hm_min == 0:
                    return np.zeros_like(hm)
                return (hm - hm_min) / (hm_max - hm_min)

            # 4. Normalize raw lenses and prepare outputs
            norm_noise = normalize_heatmap(noise_heatmap)
            norm_dct = normalize_heatmap(dct_heatmap)
            norm_freq = normalize_heatmap(freq_heatmap)
            norm_ela = normalize_heatmap(ela_heatmap)
            norm_color = normalize_heatmap(color_heatmap) 
            
            slic_overlay_rgb = mark_boundaries(img_rgb, segments, color=(1, 0, 0))
            slic_overlay_bgr = cv.cvtColor((slic_overlay_rgb * 255).astype(np.uint8), cv.COLOR_RGB2BGR)

            # Send raw arrays to UI for real-time slider manipulation
            self.finished.emit(slic_overlay_bgr, norm_noise, norm_dct, norm_freq, norm_ela, norm_color)
            
        except Exception as e:
            self.error.emit(str(e))

# --- 2. THE UI WIDGET ---
class EnhancedSplicingWidget(ToolWidget):
    def __init__(self, image, parent=None):
        super().__init__(parent)
        self.image = image
        
        # --- SMART FORENSIC DOWNSCALING ---
        max_dimension = 1080
        h, w = self.image.shape[:2]
        
        if max(h, w) > max_dimension:
            scale = max_dimension / float(max(h, w))
            new_w, new_h = int(w * scale), int(h * scale)
            # Resize the master image immediately using INTER_AREA
            self.image = cv.resize(self.image, (new_w, new_h), interpolation=cv.INTER_AREA)
        # ----------------------------------
        
        self.worker = None
        
        # Memory to store the raw math so sliders are instant
        self.raw_noise = None
        self.raw_dct = None
        self.raw_freq = None

        self.run_button = QPushButton("Compute Splicing Analysis")
        modify_font(self.run_button, bold=True)
        self.run_button.clicked.connect(self.run_analysis)

        gray_placeholder = np.full_like(self.image, 127)
        
        # --- LEFT SIDE: The Tabbed Viewers ---
        self.tabs = QTabWidget()
        self.slic_viewer = ImageViewer(self.image, gray_placeholder, "SLIC Superpixels", export=True)
        self.noise_viewer = ImageViewer(self.image, gray_placeholder, "High-Frequency Noise", export=True)
        self.dct_viewer = ImageViewer(self.image, gray_placeholder, "DCT Grid", export=True)
        self.freq_viewer = ImageViewer(self.image, gray_placeholder, "Resampling Frequency", export=True)
        self.ela_viewer = ImageViewer(self.image, gray_placeholder, "Error Level Analysis", export=True)
        self.color_viewer = ImageViewer(self.image, gray_placeholder, "Color Variance", export=True)
        
        self.tabs.addTab(self.slic_viewer, "SLIC")
        self.tabs.addTab(self.noise_viewer, "Noise")
        self.tabs.addTab(self.dct_viewer, "DCT")
        self.tabs.addTab(self.freq_viewer, "Frequency")
        self.tabs.addTab(self.ela_viewer, "ELA")
        self.tabs.addTab(self.color_viewer, "Color") 

        # --- RIGHT SIDE: The Main Heatmap ---
        self.heatmap_viewer = ImageViewer(self.image, gray_placeholder, "Final Consensus Heatmap", export=True)

        # Sync Pan/Zoom
        self.heatmap_viewer.viewChanged.connect(self.slic_viewer.changeView)
        self.heatmap_viewer.viewChanged.connect(self.noise_viewer.changeView)
        self.heatmap_viewer.viewChanged.connect(self.dct_viewer.changeView)
        self.heatmap_viewer.viewChanged.connect(self.freq_viewer.changeView)
        self.heatmap_viewer.viewChanged.connect(self.ela_viewer.changeView)
        self.heatmap_viewer.viewChanged.connect(self.color_viewer.changeView)
        
        self.slic_viewer.viewChanged.connect(self.heatmap_viewer.changeView)
        self.noise_viewer.viewChanged.connect(self.heatmap_viewer.changeView)
        self.dct_viewer.viewChanged.connect(self.heatmap_viewer.changeView)
        self.freq_viewer.viewChanged.connect(self.heatmap_viewer.changeView)
        self.ela_viewer.viewChanged.connect(self.heatmap_viewer.changeView)
        self.color_viewer.viewChanged.connect(self.heatmap_viewer.changeView)

        # --- SLIDER CONTROL PANEL ---
        controls_layout = QHBoxLayout()
        
        # Helper to create a slider + label
        def create_slider(name, default_val):
            layout = QVBoxLayout()
            label = QLabel(f"{name}: {default_val}%")
            label.setAlignment(Qt.AlignCenter)
            slider = QSlider(Qt.Horizontal)
            slider.setRange(0, 100)
            slider.setValue(default_val)
            slider.setEnabled(False) # Disabled until math finishes
            
            # Connect slider to label text AND the real-time fusion function
            slider.valueChanged.connect(lambda v: label.setText(f"{name}: {v}%"))
            slider.valueChanged.connect(self.update_fusion)
            
            layout.addWidget(label)
            layout.addWidget(slider)
            return slider, layout

        self.sld_noise, l_noise = create_slider("Noise Weight", 25)
        self.sld_dct, l_dct = create_slider("DCT Weight", 25)
        self.sld_freq, l_freq = create_slider("Freq Weight", 25)
        self.sld_ela, l_ela = create_slider("ELA Weight", 25)
        self.sld_color, l_color = create_slider("Color Weight", 25)
        self.sld_thresh, l_thresh = create_slider("Threshold", 60)
        
        controls_layout.addLayout(l_noise)
        controls_layout.addLayout(l_dct)
        controls_layout.addLayout(l_freq)
        controls_layout.addLayout(l_ela)
        controls_layout.addLayout(l_color)
        controls_layout.addLayout(l_thresh)

        # Build Main Layout
        main_layout = QGridLayout()
        main_layout.addWidget(self.tabs, 0, 0)
        main_layout.addWidget(self.heatmap_viewer, 0, 1)
        main_layout.addLayout(controls_layout, 1, 0, 1, 2) # Sliders span both columns
        main_layout.addWidget(self.run_button, 2, 0, 1, 2) 
        
        self.setLayout(main_layout)

    def run_analysis(self):
        self.start_time = time()
        self.run_button.setText("Computing Multi-Modal Lenses... Please wait.")
        modify_font(self.run_button, bold=False, italic=True)
        self.run_button.setEnabled(False)
        
        self.worker = ConsensusWorker(self.image)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)
        self.worker.start()

    def apply_colormap(self, matrix, colormap=cv.COLORMAP_VIRIDIS):
        uint8_img = (matrix * 255).astype(np.uint8)
        return cv.applyColorMap(uint8_img, colormap)

    def normalize_heatmap(self, hm):
        hm_min, hm_max = hm.min(), hm.max()
        if hm_max - hm_min == 0:
            return np.zeros_like(hm)
        return (hm - hm_min) / (hm_max - hm_min)

    def on_finished(self, slic_img, noise_mat, dct_mat, freq_mat, ela_mat, color_mat):
        try:
            # Store raw math in memory
            self.raw_noise = noise_mat
            self.raw_dct = dct_mat
            self.raw_freq = freq_mat
            self.raw_ela = ela_mat 
            self.raw_color = color_mat
            
            # Enable sliders
            self.sld_noise.setEnabled(True)
            self.sld_dct.setEnabled(True)
            self.sld_freq.setEnabled(True)
            self.sld_thresh.setEnabled(True)
            self.sld_ela.setEnabled(True)
            self.sld_color.setEnabled(True)

            # Update Left Tabs
            self.slic_viewer.update_processed(slic_img)
            self.noise_viewer.update_processed(self.apply_colormap(noise_mat, cv.COLORMAP_VIRIDIS))
            self.dct_viewer.update_processed(self.apply_colormap(dct_mat, cv.COLORMAP_VIRIDIS))
            self.freq_viewer.update_processed(self.apply_colormap(freq_mat, cv.COLORMAP_VIRIDIS))
            self.ela_viewer.update_processed(self.apply_colormap(ela_mat, cv.COLORMAP_VIRIDIS))
            self.color_viewer.update_processed(self.apply_colormap(color_mat, cv.COLORMAP_VIRIDIS))

            # --- THE MAGIC TRIGGER ---
            # Instead of manually drawing it, let PCA and K-Means take the wheel!
            self.auto_optimize_all()
            
            elapsed = time() - self.start_time
            self.run_button.setText(f"Analysis Complete ({elapsed:.1f} s)")
            modify_font(self.run_button, bold=True, italic=False)
            self.run_button.setEnabled(True)
            
            if hasattr(self, 'info_message'):
                self.info_message.emit("Enhanced Splicing Analysis Complete.")
                
        except Exception as e:
            import traceback
            QMessageBox.critical(self, "Hidden Error Caught!", f"The UI crashed at:\n{traceback.format_exc()}")
            self.run_button.setEnabled(True)

    # --- THE REAL-TIME FUSION ENGINE ---
    # This runs instantly every time you drag a slider!
    def update_fusion(self):
        if self.raw_noise is None: return # Prevent running before math is done
        
        try:
            # 1. Grab values from sliders
            w_n = self.sld_noise.value() / 100.0
            w_d = self.sld_dct.value() / 100.0
            w_f = self.sld_freq.value() / 100.0
            thresh = self.sld_thresh.value() / 100.0
            w_e = self.sld_ela.value() / 100.0
            w_c = self.sld_color.value() / 100.0

            # 2. Fuse the stored arrays together
            consensus = (self.raw_noise * w_n) + (self.raw_dct * w_d) + (self.raw_freq * w_f) + (self.raw_ela * w_e) + (self.raw_color * w_c)
            final_consensus = self.normalize_heatmap(consensus)
            
            # --- THE LUMINANCE MASK ---
            # 1. Convert original image to grayscale to check brightness
            img_gray = cv.cvtColor(self.image, cv.COLOR_BGR2GRAY)
            # 2. Create a mask where pixels darker than 240 are 1.0, and blown-out bright pixels are 0.0
            luminance_mask = (img_gray < 240).astype(np.float32)
            # 3 . Multiply the consensus map by the mask to instantly erase false bright spots!
            final_consensus = final_consensus * luminance_mask
            # -------------------------------
            
            # 3. Apply the threshold from the slider
            binary_map = (final_consensus > thresh).astype(np.uint8) * 255
            
            # 4. Morphological Cleanup
            kernel = np.ones((5, 5), np.uint8)
            cleaned_map = cv.morphologyEx(binary_map, cv.MORPH_OPEN, kernel)
            cleaned_map = cv.dilate(cleaned_map, kernel, iterations=2)
            
            # 5. Draw the Base Heatmap
            heatmap_colored = self.apply_colormap(final_consensus, cv.COLORMAP_JET)

            # 6. Find ALL contours, and apply the SMART SHAPE FILTER!
            contours, _ = cv.findContours(cleaned_map, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE)
            for c in contours:
                area = cv.contourArea(c)
                if area > 150: # Ignore tiny 10x10 pixel dust specks
                    x, y, w, h = cv.boundingRect(c)
                    
                    # --- THE SMART SHAPE FILTER ---
                    # 1. Aspect Ratio: Is it extremely long and thin? (Like a blade of grass)
                    aspect_ratio = float(w) / float(h)
                    
                    # 2. Extent: Is it a hollow, stringy shape? (Area compared to its bounding box)
                    rect_area = w * h
                    extent = float(area) / float(rect_area)
                    
                    # Rule: If it is super stringy (ratio > 4 or < 0.25) OR mostly empty space (extent < 0.25)...
                    if (aspect_ratio > 4.0 or aspect_ratio < 0.25) or (extent < 0.25):
                        continue # ...IGNORE IT! (It's probably grass or background noise)
                    # ------------------------------
                    
                    # If it survives the filter, draw the box on heatmap_colored!
                    cv.rectangle(heatmap_colored, (int(x), int(y)), (int(x + w), int(y + h)), (0, 255, 0), 3) 
                    cv.putText(heatmap_colored, "SPLICING DETECTED", (int(x), int(y) - 10), 
                               cv.FONT_HERSHEY_SIMPLEX, 0.7, (0, 255, 0), 2)
            
            self.heatmap_viewer.update_processed(heatmap_colored)
            
        except Exception as e:
            print(f"Slider Update Error: {e}")
    
    def on_error(self, err):
        QMessageBox.critical(self, "Error", f"Analysis failed: {err}")
        self.run_button.setText("Compute Splicing Analysis")
        modify_font(self.run_button, bold=True, italic=False)
        self.run_button.setEnabled(True)
        
    def auto_optimize_all(self):
        """Unsupervised ML: Uses PCA for Weights and K-Means for the Threshold"""
        if self.raw_noise is None: return

        try:
            # --- PART 1: PCA FOR DYNAMIC LENS WEIGHTING ---
            # 1. Flatten all 5 heatmaps into 1D arrays
            n, d, f = self.raw_noise.flatten(), self.raw_dct.flatten(), self.raw_freq.flatten()
            e, c = self.raw_ela.flatten(), self.raw_color.flatten()
            
            # 2. Stack them into a matrix and Standardize (Crucial for PCA)
            X = np.column_stack([n, d, f, e, c])
            X_norm = (X - np.mean(X, axis=0)) / (np.std(X, axis=0) + 1e-8)
            
            # 3. Calculate Covariance and Eigenvectors (The PCA Math)
            cov_mat = np.cov(X_norm, rowvar=False)
            _, eigen_vecs = np.linalg.eigh(cov_mat)
            
            # 4. Extract Principal Component 1 (The axis of maximum variance)
            pc1 = eigen_vecs[:, -1] 
            
            # 5. Convert PC1 into percentages for our 5 sliders
            weights_pct = (np.abs(pc1) / np.sum(np.abs(pc1))) * 100.0

            # 6. Snap the UI sliders! (Block signals to prevent UI stuttering)
            sliders = [self.sld_noise, self.sld_dct, self.sld_freq, self.sld_ela, self.sld_color]
            for sld in sliders: sld.blockSignals(True)

            self.sld_noise.setValue(int(weights_pct[0]))
            self.sld_dct.setValue(int(weights_pct[1]))
            self.sld_freq.setValue(int(weights_pct[2]))
            self.sld_ela.setValue(int(weights_pct[3]))
            self.sld_color.setValue(int(weights_pct[4]))

            for sld in sliders:
                sld.blockSignals(False)
                sld.valueChanged.emit(sld.value()) # Update the UI text labels

            # --- PART 2: BUILD CONSENSUS WITH NEW AI WEIGHTS ---
            w_n, w_d, w_f, w_e, w_c = weights_pct / 100.0
            consensus = (self.raw_noise * w_n) + (self.raw_dct * w_d) + (self.raw_freq * w_f) + (self.raw_ela * w_e) + (self.raw_color * w_c)
            final_consensus = self.normalize_heatmap(consensus)

            # --- PART 3: K-MEANS FOR AUTO-THRESHOLD ---
            pixel_values = final_consensus.reshape((-1, 1)).astype(np.float32)
            criteria = (cv.TERM_CRITERIA_EPS + cv.TERM_CRITERIA_MAX_ITER, 10, 1.0)
            _, _, centers = cv.kmeans(pixel_values, 2, None, criteria, 10, cv.KMEANS_PP_CENTERS)
            
            optimal_thresh = float((centers[0][0] + centers[1][0]) / 2.0)
            
            # 4. Snap the threshold slider!
            # Note: setting this automatically triggers update_fusion() to draw the final box!
            self.sld_thresh.setValue(int(optimal_thresh * 100))

        except Exception as err:
            print(f"Auto-Optimization Error: {err}")