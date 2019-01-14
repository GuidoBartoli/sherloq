#ifndef QUALITYWIDGET_H
#define QUALITYWIDGET_H

#include "toolwidget.h"
#include <opencv2/opencv.hpp>

class QualityWidget : public ToolWidget
{
    Q_OBJECT

public:
    QualityWidget(const QString &fileName, QWidget* parent = Q_NULLPTR);

private:
    void showError(int type);
    bool findNextMarker(FILE* file, uchar byte1, uchar byte2, uchar byte3);
    void getTables(int quality, cv::Mat &luma_qt, cv::Mat &chroma_qt, bool baseline = true);

};

#endif // QUALITYWIDGET_H
