#ifndef ECHOWIDGET_H
#define ECHOWIDGET_H

#include <opencv2/opencv.hpp>
#include "toolwidget.h"
#include "singleviewer.h"
#include "clickableslider.h"

class EchoWidget : public ToolWidget
{
    Q_OBJECT

public:
    EchoWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    QSpinBox*        radiusSpin;
    ClickableSlider* contrastSlider;
    QLabel*          contrastLabel;
    SingleViewer*    echoViewer;
    std::vector<cv::Mat> channels;
    cv::Mat computeEcho();

private slots:
    void updateView();

};

#endif // ECHOWIDGET_H
