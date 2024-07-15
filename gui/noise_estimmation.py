""" 
This code implements local noise estimation based on on high pass wavelet coeffiecients and subsequent grid blocking (median based).
The technique is explain in the following paper:
"Using noise inconsistencies for blind image forensics" by Babak Mahdian & Stanislav Saic
A block merging step has not been included because all attempts yielded unsastifactory results and made analysis more difficult.

In the paper the merging step appears highly effective, but this could be because of the isolated test conditions where the only obserable
noise in the image was the one added by the researchers for testing puposes.
"""

from PySide6.QtWidgets import QVBoxLayout, QHBoxLayout, QLabel, QSpinBox, QPushButton
from tools import ToolWidget
from viewer import ImageViewer

#algorithm nessessary imports
import numpy as np
import pywt
import cv2

class NoiseWaveletBlockingWidget(ToolWidget):
    #tool layout
    def __init__(self, filename, image, parent=None):
        super(NoiseWaveletBlockingWidget, self).__init__(parent)

        #save variables to self
        self.filename = filename
        self.image = image

        #prepare user interface - input variables
        self.blocksize_spin = QSpinBox()
        self.blocksize_spin.setRange(1, 100)
        self.blocksize_spin.setValue(8)

        self.process_button = QPushButton(self.tr("Process image"))

        #combine top layout
        top_layout = QHBoxLayout()
        top_layout.addWidget(QLabel(self.tr("Block size:")))
        top_layout.addWidget(self.blocksize_spin)
        top_layout.addWidget(self.process_button)
        top_layout.addStretch()

        self.viewer = ImageViewer(image, image, None)

        self.calculate_noise_map()
        self.process_button.clicked.connect(self.calculate_noise_map)

        #main layout
        main_layout = QVBoxLayout()
        main_layout.addLayout(top_layout)
        main_layout.addWidget(self.viewer)
        self.setLayout(main_layout)




    def calculate_noise_map(self):
        self.process_button.setEnabled(False) #wait for processing
        blocksize = self.blocksize_spin.value()

        im = cv2.imread(self.filename, cv2.IMREAD_GRAYSCALE)

        y = np.double(im)
        #3.1 wavelet transform
        cA1, (cH, cV, cD) = pywt.dwt2(y, 'db8')
        
        cD = cD[:cD.shape[0] // blocksize * blocksize, :cD.shape[1] // blocksize * blocksize]
        
        #3.2 non overlapping blocks
        block = np.zeros((cD.shape[0] // blocksize, cD.shape[1] // blocksize, blocksize ** 2))
        
        for ii in range(0, cD.shape[0] - blocksize + 1, blocksize):
            for jj in range(0, cD.shape[1] - blocksize + 1, blocksize):
                block_elements = cD[ii:ii+blocksize, jj:jj+blocksize]
                block[ii // blocksize, jj // blocksize, :] = block_elements.flatten()

        #3.3 noise level estimation
        noise_map = np.median(np.abs(block), axis=2) / 0.6745
        #3.4 blocks merging - 
        #not included, merging results for real images were dissatisfactory, merging works better in lab conditions

        #convert to viewer expectations and update viewer
        noise_map_8u = cv2.normalize(noise_map, None, 0, 255, cv2.NORM_MINMAX, dtype=cv2.CV_8U)
        resized_noise_map = cv2.resize(noise_map_8u, (self.image.shape[1], self.image.shape[0]), interpolation=cv2.INTER_NEAREST )
        noise_map_BGR = cv2.cvtColor(resized_noise_map, cv2.COLOR_GRAY2BGR)
        
        self.viewer.update_processed(noise_map_BGR)

        self.process_button.setEnabled(True) #allow new process to start