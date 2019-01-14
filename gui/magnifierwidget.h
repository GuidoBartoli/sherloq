#ifndef MAGNIFIERWIDGET_H
#define MAGNIFIERWIDGET_H

#include "toolwidget.h"
#include "singleviewer.h"
#include <opencv2/opencv.hpp>

class MagnifierWidget : public ToolWidget
{
    Q_OBJECT

    const float PERCENTILE = 0.05;

public:
    MagnifierWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    SingleViewer* magnifierViewer;
    QRadioButton* equalizeRadio;
    QRadioButton* contrast1Radio;
    QRadioButton* contrast2Radio;
    QRadioButton* retinexRadio;
    QRadioButton* lastRadio;
    cv::Mat input;
    cv::Mat autoLUT(const cv::Mat &channel, float percentile);

private slots:
    void updateView(QRect rect);
    void updateRect();
};

#endif // MAGNIFIERWIDGET_H
