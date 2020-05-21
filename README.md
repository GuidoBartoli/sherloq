<p align="center">
  <img src="logo/sherloq.png" width="600px" alt="Sherloq" />
  <br><b>An open source image forensic toolset</b>
</p>

# Introduction
"*Forensic Image Analysis is the application of image science and domain expertise to interpret the content of an image and/or the image itself in legal matters. Major subdisciplines of Forensic Image Analysis with law enforcement applications include: Photogrammetry, Photographic Comparison, Content Analysis, and Image Authentication.*" (Scientific Working Group on Imaging Technologies)

**Sherloq** is a personal research project about implementing a fully integrated environment for digital image forensics. It is not meant as an automatic tool that decide if an image is forged or not (that tool probably will never exist...), but as a companion in experimenting with various algorithms to discover potential image inconsistencies.

While many commercial solutions have unaffordable prices and are reserved to law enforcement and government agencies only, this toolset aims to be a powerful and extensible framework as a starting point for anyone interested in state-of-the-art forensic algorithms.

I strongly believe that *security-by-obscurity* is the wrong way to offer any kind of forensic service (i.e. "Using this proprietary software I guarantee you that this photo *is* pristine... and you have to trust me!"). Following the open-source philosophy, everyone should instead be able to try various techniques on their own, gain knowledge and share it to the community... even better if they contribute with code improvements! :)

# History
The first version was written in 2015 using C++11 to build a command line utility with many options, but soon it turned to be too cumbersome and not much interactive. That version is contained in the "old/cli" folder and can be compiled with CMake after installing OpenCV, Boost and AlgLib libraries. This first proof of concept offered about 80% of planned features (see below for the full list).

While also including novel algorithms, the 2017 version mainly added a Qt-based multi-window GUI to provide a better user experience. Multiple analyses could be shown on screen and a fast zoom/scroll  viewer was implemented for easier image navigation. This project is contained in the "old/gui" folder and can be compiled with Qt Creator with Qt 5 and OpenCV 3 and covered about 70% of planned features (see below for the full list).

Fast forward to 2020 when I decided to port everything in Python + PySide2 + OpenCV for a much easier installation and maintenance. This iteration is just begun and I have ported about 40% of the previous code on the new platform, but I think this will be the final "form" of the project (as long as someone does not volunteer to develop a nice web application!).

I'm happy to share my code and get in contact with anyone interested to improve or test it, but please keep in mind that this repository is *not* intended for distributing a final product, my aim is just to publicly track development of an unpretentious educational tool, so expect bugs, unpolished code and missing features! ;)

# Install
To run the Python version, just create a Python 3 virtual environment and execute `pip install -r requirements.txt` from inside the `python` folder, then `python sherloq.py`.

These are my suggested steps to create a virtual environment on Ubuntu-based distributions:

```
Install package manager
$ sudo apt install python3-distutils python3-dev subversion
$ wget https://bootstrap.pypa.io/get-pip.py
$ sudo python3 get-pip.py

Enable virtual environments
$ sudo pip install virtualenv virtualenvwrapper
$ echo -e "\n# Python Virtual Environments" >> ~/.bashrc
$ echo "export WORKON_HOME=$HOME/.virtualenvs" >> ~/.bashrc
$ echo "export VIRTUALENVWRAPPER_PYTHON=/usr/bin/python3" >> ~/.bashrc
$ echo "source /usr/local/bin/virtualenvwrapper.sh" >> ~/.bashrc
$ source ~/.bashrc

Create virtual environment
$ mkvirtualenv dev -p python3
```

Compiling the `cli` or `gui` version can be *tricky*, since it is old code and tested only on Ubuntu-based distributions. When everything will be functional in Python it will be very likely removed, however, if you need assistance to set it up, drop me a line and I will try to look into itm (no guarantees!).

# Features
This list contains the functions that the toolkit could provide once the beta stage is reached.

## General
- __Original Image__: display the unaltered reference image for visual inspection
- __File Digest__: retrieve file information and compute many hashes and ballistics
- __Similarity Search__: use reverse search services for finding similar images online
- __Automatic Tagging__: apply deep learning algorithms for automatic picture tagging

## Metadata
- __Header Structure__: dump the physical EXIF structure and display an interactive view
- __Metadata Extraction__: scan through file metadata and gather all available information
- __Thumbnail Analysis__: extract optional embedded thumbnail and compare with original
- __Geolocation Data__: retrieve optional geo-location data and show it on a world map

## Inspection
- __Enhancing Magnifier__: use various visual enhancement for better identifying forgeries
- __Reference Comparison__: open a synchronized double view to compare two different pictures
- __Global Adjustments__: apply standard adjustments (contrast, brightness, hue, saturation)
- __Range Compression__: compress tonality over different ranges to detect inconsistencies

## JPEG
- __Quality Estimation__: extract quantization tables and estimate last saved JPEG quality
- __Error Level Analysis__: show pixel-level difference against different compression levels
- __Quantization Ghosts__: Use residuals to detect multiple compressions at different levels
- __Double Compression__: exploit DCT First-Digit-Statistics to detect double compression

## Colors
- __RGB/HSV 3D Plots__: display interactive 2D and 3D plots of RGB and HSV pixel data
- __RGB PCA Projection__: use color PCA to project RGB values onto reduced vector spaces
- __RGB Pixel Statistics__: compute minimum/maximum/average RGB values for every pixel
- __Color Space Conversion__: convert color channels into RGB/HSV/YCbCr/Lab/CMYK color spaces

## Tonality
- __Luminance Gradient__: analyze horizontal and vertical brightness variations of the image
- __Echo Edge Filter__: use derivative filter to reveal artificial out-of-focus zones
- __Frequency Separation__: estimate high and low frequency components of the luminance channel
- __Wavelet Reconstruction__: reconstruct image with different wavelet coefficient thresholds

## Noise
- __Noise Extraction__: estimate and visualize gaussian noise components of the image
- __Min/Max Deviation__: highlight pixels deviating from block-based min/max statistics
- __SNR Consistency__: evaluate uniformity of signal-to-noise ratio across the image
- __Noise Segmentation__: cluster noise into uniform regions for anomaly detection

## Tampering
- __Contrast Enhancement__: analyze color distribuions to detect contrast enhancements
- __Region Cloning__: use feature descriptors for copy/rotate clone area detection
- __Image Resampling__: analyze 2D pixel interpolation for detecting resampling traces
- __Composite Splicing__: exploit DCT statistics for automatic splicing zone detection

# Screenshots
Here are some screenshots from the previous C++ Qt GUI:
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

# Bibliography
- Black Hat Briefings DC. (2008) A Picture's Worth: Digital Image Analysis and Forensics [White paper]. Washington, DC. Retrieved from http://blackhat.com/presentations/bh-dc-08/Krawetz/Whitepaper/bh-dc-08-krawetz-WP.pdf
