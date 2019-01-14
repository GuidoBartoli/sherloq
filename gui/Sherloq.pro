#-------------------------------------------------
#
# Project created by QtCreator 2017-04-28T13:45:32
#
#-------------------------------------------------

QT += core gui webenginewidgets charts
#qtHaveModule(opengl): QT += opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET    = Sherloq
TEMPLATE  = app
CONFIG   += c++11

SOURCES  += main.cpp\
    mainwindow.cpp \
    dynamicview.cpp \
    digestwidget.cpp \
    ExifTool.cpp \
    ExifToolPipe.cpp \
    TagInfo.cpp \
    metadatawidget.cpp \
    structurewidget.cpp \
    elawidget.cpp \
    utility.cpp \
    cvmatandqimage.cpp \
    singleviewer.cpp \
    toolwidget.cpp \
    originalwidget.cpp \
    tooltree.cpp \
    thumbnailwidget.cpp \
    locationwidget.cpp \
    qualitywidget.cpp \
    ghostwidget.cpp \
    gradientwidget.cpp \
    sweepwidget.cpp \
    rgbpixelwidget.cpp \
    minmaxwidget.cpp \
    spacewidget.cpp \
    colorpcawidget.cpp \
    adjustmentwidget.cpp \
    clickableslider.cpp \
    separationwidget.cpp \
    magnifierwidget.cpp \
    doublewidget.cpp \
    comparisonwidget.cpp \
    echowidget.cpp \
    contrastwidget.cpp \
    alglib/alglibinternal.cpp \
    alglib/alglibmisc.cpp \
    alglib/ap.cpp \
    alglib/dataanalysis.cpp \
    alglib/diffequations.cpp \
    alglib/fasttransforms.cpp \
    alglib/integration.cpp \
    alglib/interpolation.cpp \
    alglib/linalg.cpp \
    alglib/optimization.cpp \
    alglib/solvers.cpp \
    alglib/specialfunctions.cpp \
    alglib/statistics.cpp \
    snrwidget.cpp \
    histogramwidget.cpp

HEADERS  += mainwindow.h \
    dynamicview.h \
    digestwidget.h \
    ExifTool.h \
    ExifToolPipe.h \
    TagInfo.h \
    metadatawidget.h \
    structurewidget.h \
    elawidget.h \
    utility.h \
    cvmatandqimage.h \
    singleviewer.h \
    toolwidget.h \
    originalwidget.h \
    tooltree.h \
    thumbnailwidget.h \
    locationwidget.h \
    qualitywidget.h \
    ghostwidget.h \
    gradientwidget.h \
    sweepwidget.h \
    rgbpixelwidget.h \
    minmaxwidget.h \
    spacewidget.h \
    colorpcawidget.h \
    adjustmentwidget.h \
    clickableslider.h \
    separationwidget.h \
    magnifierwidget.h \
    doublewidget.h \
    comparisonwidget.h \
    echowidget.h \
    contrastwidget.h \
    alglib/alglibinternal.h \
    alglib/alglibmisc.h \
    alglib/ap.h \
    alglib/dataanalysis.h \
    alglib/diffequations.h \
    alglib/fasttransforms.h \
    alglib/integration.h \
    alglib/interpolation.h \
    alglib/linalg.h \
    alglib/optimization.h \
    alglib/solvers.h \
    alglib/specialfunctions.h \
    alglib/statistics.h \
    alglib/stdafx.h \
    snrwidget.h \
    histogramwidget.h

#INCLUDEPATH += /opt/ros/kinetic/include/opencv-3.2.0-dev/
#DEPENDPATH  += /opt/ros/kinetic/include/opencv-3.2.0-dev/
#LIBS        += \
#    -L/opt/ros/kinetic/lib/ \
#    -lopencv_core3      \
#    -lopencv_imgproc3   \
#    -lopencv_imgcodecs3 \
#    -lopencv_ml3

LIBS     += \
    -L/usr/local/lib/ \
    -lopencv_core    \
    -lopencv_imgproc \
    -lopencv_imgcodecs
    -lopencv_ml
