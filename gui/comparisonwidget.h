#ifndef COMPARISONWIDGET_H
#define COMPARISONWIDGET_H

#include "toolwidget.h"
#include "singleviewer.h"
#include <opencv2/opencv.hpp>

class ComparisonWidget : public ToolWidget
{
    Q_OBJECT

public:
    ComparisonWidget(const cv::Mat &evidence, QWidget* parent = Q_NULLPTR);

private:
    cv::Mat evidence;
    cv::Mat reference;
    cv::Mat difference;
    cv::Mat indexMap;
    SingleViewer* referenceViewer;
    SingleViewer* evidenceViewer;
    QRadioButton* originalRadio;
    QRadioButton* differenceRadio;
    QRadioButton* indexMapRadio;
    QLabel* filenameLabel;
    QLabel* indicesLabel;

    double ssim(const cv::Mat &img_src, const cv::Mat &img_compressed, int block_size = 8);
    double sigma(const cv::Mat &m, int i, int j, int block_size);
    double cov(const cv::Mat &m1, const cv::Mat &m2, int i, int j, int block_size);
    double eqm(const cv::Mat &img1, const cv::Mat &img2);
    double psnr(const cv::Mat &img1, const cv::Mat &img2);

    double ssim2(const cv::Mat &img1orig, const cv::Mat &img2orig);
private slots:
    void loadReference();\
    void updateReference();

};

#endif // COMPARISONWIDGET_H
