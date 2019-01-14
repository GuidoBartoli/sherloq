#ifndef CONTRASTWIDGET_H
#define CONTRASTWIDGET_H

#include "toolwidget.h"
#include <opencv2/opencv.hpp>

class ContrastWidget : public ToolWidget
{
    Q_OBJECT

public:
    ContrastWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    double  computeHistErr(const cv::Mat &input);
    double  computeChanSim(const cv::Mat &input);
    cv::Mat histDFT(const cv::Mat &hist, bool full = false);
    void contrastEnhancement(const cv::Mat &input, cv::Mat &output, QWidget *parent);
    void contrastEnhancement2(const cv::Mat &input, cv::Mat &output, QWidget *parent);
};

#endif // CONTRASTWIDGET_H
