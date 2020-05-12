<p align="center">
  <img src="logo/sherloq.png" width="600px" alt="Sherloq" />
  <br><b>An open source image forensic toolset</b>
</p>

# Introduction
"*Forensic Image Analysis is the application of image science and domain expertise to interpret the content of an image and/or the image itself in legal matters. Major subdisciplines of Forensic Image Analysis with law enforcement applications include: Photogrammetry, Photographic Comparison, Content Analysis, and Image Authentication.*" (Scientific Working Group on Imaging Technologies)

**Sherloq** is a personal research project about implementing a fully integrated environment for digital image forensics. It is not meant as an automatic tool that decide if an image is forged or not (that tool probably will never exist...), but as a companion in experimenting with various algorithms to discover potential image inconsistencies.

While many commercial solutions have unaffordable prices and are reserved to law enforcement and government agencies only, this toolset aims to be a powerful and extensible framework as a starting point for anyone interested in state-of-the-art forensic algorithms.

I strongly believe that *security-by-obscurity* is the wrong way to offer any kind of forensic service (i.e. "Using this proprietary software I guarantee you that this photo *is* pristine... and you have to trust me!"). Instead, following the open-source philosophy, everyone should be able to try various techniques on their own, gain knowledge and share it to the community... even better if they contribute with code improvements! :)

# History
The first version was written in 2015 using C++11 to build a command line utility with many options, but soon it turned to be too complicated and not much interactive. That version is contained in the "cli" folder and can be compiled with CMake after installing OpenCV, Boost and AlgLib libraries. This first proof of concept offered about 80% of planned features (see below).

While also including novel algorithms, the 2017 version mainly added a Qt-based multi-window GUI to provide a better user experience. Multiple analyses could be shown on screen and a fast zoom/scroll  viewer was implemented for easier image navigation. This project is contained in the "gui" folder and can be compiled with Qt Creator with Qt 5 and OpenCV 3 and covered about 75% of planned features (see below).

Fast forward to 2020 when I decided to port everything in Python + PySide2 + OpenCV for a much easier installation and maintenance. This iteration is just begun and I have ported about 15% of the previous code on the new platform, but I think this will be the final "form" of the project (as long as someone does not volunteer to develop a web application!)

I'm happy to share my code and get in contact with anyone interested to improve or test it, but please keep in mind that this repository is *not* intended for distributing a final product, my intent is just to publicly track development, so expect bugs, unpolished code and missing features! ;) 

# Install
Compiling the `cli` or `gui` version can be tricky, since it is old code and tested only on Ubuntu-based distributions, however, if you need assistance to set it up, drop me a line and I will try to look into it (no guarantees!).

To run the Python version, just create a virtual environment and execute `pip install -r requirements.txt` from inside the `python` folder, then `python sherloq.py`.

# Features
This list contains all functions that **Sherloq** could provide once the beta stage is reached.

## General
- __Original Image__: display the unaltered reference image for visual inspection (&ast;&ast;&ast;)
- __Image Digest__: compute various hashes together with extension ballistics (&ast;&ast;)
- __Similarity Search__: use reverse search services for finding similar images on the web (&ast;)
- __Automatic Tagging__: apply deep learning algorithms for automatic picture tagging (&ast;)

## Metadata
- __Header Structure__: dump the physical EXIF structure and display an interactive view (&ast;&ast;&ast;)
- __Metadata Extraction__: gather all metadata information and display security warnings (&ast;&ast;)
- __Thumbnail Analysis__: if present, extract embedded thumbnail and highlight discrepancies (&ast;&ast;&ast;)
- __Geolocation Data__: if present, get geographic data and locate them on a world map view (&ast;&ast;&ast;)

## Inspection
- __Enhancing Magnifier__: apply local visual enhancements for better identifying forgeries (&ast;&ast;&ast;)
- __Image Adjustments__: apply standard adjustments (contrast, brightness, hue, saturation, ...) (&ast;&ast;&ast;)
- __Tonal Range Sweep__: interactive tonality range compression for easier artifact detection (&ast;&ast;&ast;)
- __Reference Comparison__: synchronized double view to compare reference and evidence images (&ast;&ast;&ast;)

## JPEG
- __Quality Estimation__: extract quantization tables and estimate last saved JPEG quality (&ast;&ast;&ast;)
- __Compression Ghosts__: use error residuals to detect multiple compressions at different levels (&ast;&ast;)
- __Double Compression__: exploit First Digit Statistics to discover potential double compression (&ast;&ast;)
- __Error Level Analysis__: identify areas with different compression levels against a fixed quality (&ast;&ast;&ast;)

## Colors
- __RGB/HSV 3D Plots__: display interactive 2D and 3D plots of RGB and HSV pixel data (&ast;)
- __Color Space Conversion__: convert image into RGB/HSV/YCbCr/Lab/CMYK color spaces (&ast;&ast;&ast;)
- __Principal Component Analysis__: use PCA to project RGB values onto a different vector space (&ast;&ast;&ast;)
- __RGB Pixel Statistics__: compute minimum/maximum/average RGB values for every pixel (&ast;&ast;&ast;)

## Tonality
- __Luminance Gradient__: analyze brightness variations along X/Y axes of the image (&ast;&ast;&ast;)
- __Frequency Separation__: extract the finest details of the luminance channel (&ast;)
- __Echo Edge Filter__: use 2D Laplacian filter to reveal artificial blurred zones (&ast;&ast;&ast;)
- __Wavelet Reconstruction__: re-synthesize image varying wavelet coefficient thresholds (&ast;)

## Noise
- __Noise Extraction__: estimate and separate the natural noise component of the image (&ast;&ast;&ast;)
- __Min/Max Deviation__: highlight pixels deviating from block-based min/max statistics (&ast;&ast;&ast;)
- __SNR Consistency__: evaluate uniformity of signal-to-noise ratio across the image (&ast;&ast;&ast;)
- __Noise Segmentation__: cluster uniform noise areas for easier tampering detection (&ast;)

## Tampering
- __Contrast Enhancement__: analyze histogram inconsistencies caused by enhancements (&ast;&ast;&ast;)
- __Clone Detection__: use invariant feature descriptors for copy/rotate clone area detection (&ast;&ast;)
- __Resampling Detection__: analyze 2D pixel interpolation for detecting resampling traces (&ast;&ast;)
- __Splicing Detection__: use DCT coefficient statistics for automatic splicing zone detection (&ast;) 

# Screenshots
<p align="center">
  <img src="screenshots/File.jpg" alt="File Analysis"/>
  <br><b>File Analysis</b>: Metadata, Digest and EXIF
</p>

<p align="center">
  <img src="screenshots/Color.jpg" alt="Color Analysis"/>
  <br><b>Color Analysis</b>: Space Conversion, PCA Projection, Histograms and Statistics
</p>

<p align="center">
  <img src="screenshots/Visual.jpg" alt="Visual Inspection"/>
  <br><b>Visual Inspection</b>: Magnifier Loupe, Image Adjustments and Evidence Comparison
</p>

<p align="center">
  <img src="screenshots/JPEG.jpg" alt="JPEG Analysis"/>
  <br><b>JPEG Analysis</b>: Quantization Tables, Compression Ghosts and Error Level Analysis
</p>

<p align="center">
  <img src="screenshots/LumaNoise.jpg" alt="Luminance/Noise"/>
  <br><b>Luminance and Noise</b>: Light Gradient, Echo Edge, Min/Max Deviation and SNR Consistency
</p>
