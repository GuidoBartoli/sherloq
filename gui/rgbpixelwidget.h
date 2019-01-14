#ifndef RGBPIXELWIDGET_H
#define RGBPIXELWIDGET_H

#include "toolwidget.h"
#include "singleviewer.h"
#include <opencv2/opencv.hpp>

class RgbPixelWidget : public ToolWidget
{
    Q_OBJECT
public:
    RgbPixelWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    QRadioButton* minButton;
    QRadioButton* maxButton;
    QRadioButton* avgButton;
    QRadioButton* lastRadio = Q_NULLPTR;
    SingleViewer* rgbPixelViewer;
    cv::Mat minImage, maxImage, avgImage;

private slots:
    void updateView(bool);
};

#endif // RGBPIXELWIDGET_H
