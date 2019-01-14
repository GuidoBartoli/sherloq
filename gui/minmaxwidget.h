#ifndef MINMAXWIDGET_H
#define MINMAXWIDGET_H

#include "toolwidget.h"
#include <opencv2/opencv.hpp>

class MinMaxWidget : public ToolWidget
{
    Q_OBJECT

public:
    MinMaxWidget(const cv::Mat &input, QWidget *parent = 0);

};

#endif // MINMAXWIDGET_H
