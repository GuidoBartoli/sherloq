#ifndef ADJUSTMENTWIDGET_H
#define ADJUSTMENTWIDGET_H

#include "toolwidget.h"
#include "singleviewer.h"
#include <opencv2/opencv.hpp>
#include "clickableslider.h"

class AdjustmentWidget : public ToolWidget
{
    Q_OBJECT

public:
    AdjustmentWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    cv::Mat          input, hsv;
    SingleViewer*    adjustmentViewer;
    ClickableSlider* brightnessSlider;
    QLabel*          brightnessLabel;
    ClickableSlider* contrastSlider;
    QLabel*          contrastLabel;
    ClickableSlider* hueSlider;
    QLabel*          hueLabel;
    ClickableSlider* saturationSlider;
    QLabel*          saturationLabel;
    QCheckBox*       invertCheck;
    QPushButton*     resetButton;

private slots:
    void resetSliders();
    void updateView();
    void resetBrightness();
    void resetContrast();
    void resetHue();
    void resetSaturation();

};

#endif // ADJUSTMENTWIDGET_H
