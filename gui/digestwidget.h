#ifndef DIGESTWIDGET_H
#define DIGESTWIDGET_H

#include "toolwidget.h"
#include <opencv2/opencv.hpp>

class DigestWidget : public ToolWidget
{
    Q_OBJECT

public:
    DigestWidget(const QString &fileName, QWidget* parent = Q_NULLPTR);

private slots:
    void copyCell(QTableWidgetItem* item);

private:
    void    computePHash(const QString &fileName, QString &bytes, QString &bits);
    QString nameBallistics(const QString &fileName);
    //double  medianMat(cv::Mat Input);

};

#endif // DIGESTWIDGET_H
