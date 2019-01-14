#ifndef DYNAMICVIEW_H
#define DYNAMICVIEW_H

#include <QtWidgets>

class DynamicView : public QGraphicsView
{
    Q_OBJECT

public:
    DynamicView(QGraphicsScene* scene, QWidget* parent = Q_NULLPTR);
    //double getZoomLevel() const;
    //QRect  getViewRect() const;
    //void   setZoomFactor(double factor);
    //double getZoomFactor() const;
    void  changeView(QRect, double zoomFactor, int horizScroll, int vertScroll);
    QRect getSceneRect();

protected:
    void wheelEvent(QWheelEvent* event)        Q_DECL_OVERRIDE;
    void mouseDoubleClickEvent(QMouseEvent*)   Q_DECL_OVERRIDE;
    void mouseMoveEvent(QMouseEvent* event)    Q_DECL_OVERRIDE;
    void mousePressEvent(QMouseEvent* event)   Q_DECL_OVERRIDE;
    void mouseReleaseEvent(QMouseEvent* event) Q_DECL_OVERRIDE;
    void resizeEvent(QResizeEvent* event)      Q_DECL_OVERRIDE;

public slots:
    void fullZoom();
    void fitZoom();

private:
    void changeZoom(int direction);
    void setZoomFactor(double factor);
    void notifyViewChanged();

    bool   nextFit;
    bool   mousePressed;
    double fitFactor;

    int count = 0;

signals:
    void viewChanged(QRect sceneRect, double zoomFactor, int horizScroll, int vertScroll);
};

#endif // DYNAMICVIEW_H
