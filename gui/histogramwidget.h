#ifndef HISTOGRAMWIDGET_H
#define HISTOGRAMWIDGET_H

#include <QtWidgets>
#include <QtCharts>
#include <opencv2/opencv.hpp>
#include "toolwidget.h"

class HistogramWidget : public ToolWidget
{
    Q_OBJECT

public:
    HistogramWidget(const cv::Mat &image, QWidget* parent = 0);

private:
    QAreaSeries*   lumArea;
    QAreaSeries*   redArea;
    QAreaSeries*   greenArea;
    QAreaSeries*   blueArea;
    QChart*        histChart;

private slots:
    void setLumHist(bool toggled);
    void setRedHist(bool toggled);
    void setGreenHist(bool toggled);
    void setBlueHist(bool toggled);
    void setCompHist(bool toggled);
    void setLogHist(bool checked);

};

#endif // HISTOGRAMWIDGET_H
