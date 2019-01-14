#ifndef SINGLEVIEWER_H
#define SINGLEVIEWER_H

#include <QtWidgets>
#include <opencv2/opencv.hpp>
#include "dynamicview.h"

class SingleViewer : public QWidget
{
    Q_OBJECT

public:
    SingleViewer(const cv::Mat &original, const cv::Mat &processed = cv::Mat(),
                 const QString &title = QString(), QWidget* parent = Q_NULLPTR);
    void  updateProcessed(const cv::Mat &processed);
    void  updateOriginal(const cv::Mat &original);
    QRect getSceneRect();

signals:
    void viewChanged(QRect sceneRect, double zoomFactor, int horizScroll, int vertScroll);

public slots:
    void changeView(QRect, double zoomFactor, int horizScroll, int vertScroll);

protected:
    void keyPressEvent(QKeyEvent* event) Q_DECL_OVERRIDE;

private slots:
    void toggleMode(bool origToggled);
    void forwardViewChanged(QRect sceneRect, double zoomFactor, int horizScroll, int vertScroll);
    void exportImage();

private:
    cv::Mat original;
    cv::Mat processed;
    QLabel* zoomLabel;
    DynamicView* imageView;
    QGraphicsPixmapItem* imageItem;
    QRadioButton* originalRadio;
    QRadioButton* processRadio;

};

#endif // SINGLEVIEWER_H
