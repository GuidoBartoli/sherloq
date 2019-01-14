#ifndef DOUBLEWIDGET_H
#define DOUBLEWIDGET_H

#include "toolwidget.h"
#include <opencv2/opencv.hpp>

class DoubleWidget : public ToolWidget
{
    Q_OBJECT

public:
    DoubleWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    void firstDigitFeatures(const cv::Mat input, cv::Mat &features);
    int getFirstDigit(int number);

};

#endif // DOUBLEWIDGET_H
