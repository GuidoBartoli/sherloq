#ifndef SPACEWIDGET_H
#define SPACEWIDGET_H

#include "toolwidget.h"
#include "singleviewer.h"
#include <opencv2/opencv.hpp>

class SpaceWidget : public ToolWidget
{
    Q_OBJECT

public:
    SpaceWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    std::vector<cv::Mat> rgb;
    std::vector<cv::Mat> cmyk;
    std::vector<cv::Mat> gray;
    std::vector<cv::Mat> hsv;
    std::vector<cv::Mat> hls;
    std::vector<cv::Mat> ycrcb;
    std::vector<cv::Mat> xyz;
    std::vector<cv::Mat> lab;
    std::vector<cv::Mat> luv;

    SingleViewer* colorViewer;
    QRadioButton* rgbRadio;
    QComboBox*    rgbCombo;
    QRadioButton* cmykRadio;
    QComboBox*    cmykCombo;
    QRadioButton* grayRadio;
    QComboBox*    grayCombo;
    QRadioButton* hsvRadio;
    QComboBox*    hsvCombo;
    QRadioButton* hlsRadio;
    QComboBox*    hlsCombo;
    QRadioButton* ycrcbRadio;
    QComboBox*    ycrcbCombo;
    QRadioButton* xyzRadio;
    QComboBox*    xyzCombo;
    QRadioButton* labRadio;
    QComboBox*    labCombo;
    QRadioButton* luvRadio;
    QComboBox*    luvCombo;
    QRadioButton* lastRadio = Q_NULLPTR;


private slots:
    void updateView();
};

#endif // SPACEWIDGET_H
