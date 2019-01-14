#ifndef ORIGINALWIDGET_H
#define ORIGINALWIDGET_H

#include <QtWidgets>
#include <opencv2/opencv.hpp>
#include "toolwidget.h"

class OriginalWidget : public ToolWidget
{
    Q_OBJECT

public:
    OriginalWidget(const QString &fileName, const cv::Mat &input, QWidget* parent = Q_NULLPTR);

};

#endif // ORIGINALWIDGET_H
