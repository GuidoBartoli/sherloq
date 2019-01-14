#ifndef GRADIENTWIDGET_H
#define GRADIENTWIDGET_H

#include "toolwidget.h"
#include "singleviewer.h"
#include <opencv2/opencv.hpp>

class GradientWidget : public ToolWidget
{
    Q_OBJECT

public:
    GradientWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    cv::Mat       input;
    QSpinBox*     intensitySpin;
    QCheckBox*    invertCheck;
    QCheckBox*    modulusCheck;

    SingleViewer* gradientViewer;
    cv::Mat       scharrX, scharrY, scharrZ;
    cv::Mat       computeGradient();

private slots:
    void updateView(int);
    void resetParams();

};

#endif // GRADIENTWIDGET_H
