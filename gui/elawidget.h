#ifndef ELAWIDGET_H
#define ELAWIDGET_H

#include <opencv2/opencv.hpp>
#include "toolwidget.h"
#include "singleviewer.h"

class ElaWidget : public ToolWidget
{
    Q_OBJECT

    const int DEFAULT_QUALITY = 75;
    const int DEFAULT_SCALE   = 20;

public:
    ElaWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private slots:
    void updateView(int);
    void resetParams();

private:
    SingleViewer* elaViewer;
    QSpinBox* qualitySpin;
    QSpinBox* scaleSpin;
    cv::Mat input;
    cv::Mat computeEla();
};

#endif // ELAWIDGET_H
