#ifndef SEPARATIONWIDGET_H
#define SEPARATIONWIDGET_H

#include <opencv2/opencv.hpp>
#include "toolwidget.h"
#include "singleviewer.h"
#include "clickableslider.h"

class SeparationWidget : public ToolWidget
{
    Q_OBJECT

    const int   BLOCK_SIZE = 9;
    const int   BLOCK_PAD  = 4;
    const float MAX_STDDEV = 135;
    const int   NOISE_CONT = 80;

public:
    SeparationWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    SingleViewer*    noiseViewer;
    ClickableSlider* intensitySlider;
    QCheckBox*       accurateCheck;
    QLabel*          intensityLabel;

    cv::Mat gray, noise, noise0;
    cv::Mat computeNoise();

private slots:
    void updateView();
    void resetParams();
    void resetIntensity();

};

#endif // SEPARATIONWIDGET_H
