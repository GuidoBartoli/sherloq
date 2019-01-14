#ifndef SWEEPWIDGET_H
#define SWEEPWIDGET_H

#include "toolwidget.h"
#include "singleviewer.h"
#include <opencv2/opencv.hpp>

class SweepWidget : public ToolWidget
{
    Q_OBJECT

public:
    SweepWidget(cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    cv::Mat       input;
    QSlider*      sweepSlider;
    QSpinBox*     widthSpin;
    QLabel*       sweepLabel;
    SingleViewer* sweepViewer;
    cv::Mat       computeSweep();

private slots:
    void updateView(int);
    void resetParams();

};

#endif // SWEEPWIDGET_H
