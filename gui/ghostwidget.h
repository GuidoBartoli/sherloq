#ifndef GHOSTWIDGET_H
#define GHOSTWIDGET_H

#include "toolwidget.h"
#include <opencv2/opencv.hpp>
#include <QtCharts>

class GhostWidget : public ToolWidget
{
    Q_OBJECT

    const int DCT_SIZE   = 8;
    const int TABLE_SIZE = 64;
    const int MAX_COEFF  = 64;
    const int CHAR_SHIFT = 128;

    const int MIN_QUALITY  = 3;
    const int MAX_QUALITY  = 97;
    const int QUALITY_STEP = 1;
    const int KERNEL_SIZE  = 5;

    const float MIN_PEAK = 0.218;
    const float THR_PEAK = 0.521;
    const float MAX_PEAK = 0.823;

public:
    GhostWidget(const cv::Mat &input, QWidget* parent = Q_NULLPTR);

private:
    void computeGhosts(const cv::Mat &image, cv::Mat &qualities,
                       cv::Mat &losses, int &n, int &l1, float &p1,
                       int &l2, float &p2, int &q1, int &q2, QWidget* parent);
    void computeHist(const cv::Mat &image, cv::Mat &histogram, cv::Mat &spectrum);
    cv::Mat histogram;
    cv::Mat spectrum;
    QChart* dctChart;
    QChart* dftChart;

private slots:
    void changeMode(QTableWidgetItem* item);

};

#endif // GHOSTWIDGET_H
