#ifndef SNRWIDGET_H
#define SNRWIDGET_H

#include "toolwidget.h"
#include <opencv2/opencv.hpp>

class SnrWidget : public ToolWidget
{
    Q_OBJECT

public:
    SnrWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    void minMaxDeviation(const cv::Mat &input, cv::Mat &output);
};

#endif // SNRWIDGET_H
