""" 
    This file implements a tool for integrating TruFor into Sherloq, allowing for straightforward use of TruFor's capabilities. 
    TruFor is an AI-driven solution for digital image forensics. While powerful, AI-based approaches may not always be as reliable 
    as their statistical performance metrics might suggest. Nevertheless, TruFor can assist forensic analysts and provide 
    evidence regarding the authenticity of digital images.

    Original TruFor Project: https://github.com/grip-unina/TruFor

    original trufor work:
    https://github.com/grip-unina/TruFor

    Reference Bibtex:
    @InProceedings{Guillaro_2023_CVPR,
        author    = {Guillaro, Fabrizio and Cozzolino, Davide and Sud, Avneesh and Dufour, Nicholas and Verdoliva, Luisa},
        title     = {TruFor: Leveraging All-Round Clues for Trustworthy Image Forgery Detection and Localization},
        booktitle = {Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR)},
        month     = {June},
        year      = {2023},
        pages     = {20606-20615}
    }

    For more information about the reliability of AI in digital image forensics, we recommend reading section 3.5.5 of the following thesis:
    "Digital Image Forensics: A quantitative & qualitative comparison between State-of-the-art-AI and Traditional Techniques for detection and localization of image manipulations"
    https://github.com/UHstudent/digital_image_forensics_thesis/blob/main/Thesis%20text_Digital%20Image%20Forensics-A%20Comparative%20Study%20between%20AI%20and%20traditional%20approaches.pdf

 """

from PySide6.QtWidgets import QVBoxLayout, QHBoxLayout, QLabel
from tools import ToolWidget
from viewer import ImageViewer
import os
import sys
import matplotlib.pyplot as plt
import cv2 as cv
if os.path.exists('TruFor_main/test_docker/weights/trufor.pth.tar'):
    from TruFor_main.test_docker.src.analyze_image import process_image

class TruForWidget(ToolWidget):
    def __init__(self, filename, image, parent=None):
        super(TruForWidget, self).__init__(parent)
        # If weights exist, run TruFor process and get results
        # If the weights file does not exist, show instructions to the user
        weights_path = 'TruFor_main/test_docker/weights/trufor.pth.tar'
        if os.path.exists(weights_path):
            detection_score, localizationMap = self.trufor_process(filename)
            print(type(image))
            print(type(localizationMap))

            # Handle memory overflow or other issues
            if detection_score == -9000:
                params_layout = QHBoxLayout()
                params_layout.addWidget(QLabel(self.tr("Torch CUDA out of memory error: This error occurs when the model requires more memory than is available on your GPU. - or it could be something else, check the command prompt for the exact error message.")))
                params_layout.addStretch()
                viewer = ImageViewer(image, image)
            else:
                # Display detection information and result images
                params_layout = QHBoxLayout()
                params_layout.addWidget(QLabel(self.tr("Manipulation probability according to TruFor:")))
                params_layout.addWidget(QLabel(self.tr(str(detection_score*100) + '%')))
                params_layout.addStretch()
                viewer = ImageViewer(image, localizationMap)

            # Construct main layout
            main_layout = QVBoxLayout()
            main_layout.addLayout(params_layout)
            main_layout.addWidget(viewer)
            self.setLayout(main_layout)

        else:
            main_layout = QVBoxLayout()
            main_layout.addWidget(QLabel(self.tr("To use TruFor within Sherloq, follow these steps: "
                                                 "\n1. Download the TruFor project from https://github.com/grip-unina/TruFor."
                                                 "\n2. Rename the project folder to 'TruFor_main' (instead of 'TruFor-main')."
                                                 "\n3. Download the model weights and place them in the 'test_docker' directory."
                                                 "\nEnsure the file structure is as follows: TruFor_main\\test_docker\\weights\\trufor.pth.tar."
                                                 "\n4. Copy the 'TruFor_main' project folder into the sherloq GUI directory and re-open sherloq to enable TruFor."
                                                 "\n\nIf this message continues to appear, it likely means the file structure of the TruFor project has not been set up correctly within Sherloq.")))
            self.setLayout(main_layout)

    def trufor_process(self, filename):
        original_cwd = os.getcwd()

        try:
            # Change working directory and append paths to sys for resolving imports
            os.chdir('TruFor_main/test_docker/src/')
            sys.path.append(os.getcwd())
            sys.path.append(os.path.join(os.getcwd(), 'models'))

            gpu = 0  # GPU device, use -1 for CPU

            # Try to process the image using TruFor, handle any errors
            try:
                prediction, det_score = process_image(filename, gpu)
            except Exception as e:
                print(f"Error during process_image: {e}")
                return -9000, None  # Return error code for memory overflow or other issues

        finally:
            # Always change back to the original working directory
            os.chdir(original_cwd)

        # Save prediction temporarily and convert to BGR format for sherloq viewer
        temp_filename = "temp_mask.png"
        plt.imsave(temp_filename, prediction, cmap='RdBu_r', vmin=0, vmax=1)

        # Load the saved image as BGR using OpenCV
        prediction = cv.imread(temp_filename, cv.IMREAD_COLOR)

        # Remove the temporary file
        os.remove(temp_filename)

        # Return detection score and the localization map
        return det_score, prediction
