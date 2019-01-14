#ifndef COLORPCAWIDGET_H
#define COLORPCAWIDGET_H

#include "toolwidget.h"
#include "singleviewer.h"
#include <opencv2/opencv.hpp>

class ColorPcaWidget : public ToolWidget
{
    Q_OBJECT

public:
    ColorPcaWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    std::vector<cv::Mat> components;
    SingleViewer* pcaViewer;
    QComboBox*    componentCombo;
    QRadioButton* distanceRadio;
    QRadioButton* distanceEqRadio;
    QRadioButton* crossRadio;
    QRadioButton* lastRadio;

private slots:
    void updateView();

};

#endif // COLORPCAWIDGET_H
