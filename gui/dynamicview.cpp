#include "dynamicview.h"

//#ifndef QT_NO_OPENGL
//#include <QtOpenGL>
//#endif

DynamicView::DynamicView(QGraphicsScene* scene, QWidget* parent)
    : QGraphicsView(scene, parent)
{
    setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
    setDragMode(QGraphicsView::ScrollHandDrag);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setViewportUpdateMode(QGraphicsView::SmartViewportUpdate);
    setOptimizationFlags(QGraphicsView::DontSavePainterState);
    setRenderHint(QPainter::SmoothPixmapTransform);
    mousePressed = false;
    fitZoom();

//    connect(horizontalScrollBar(), &QScrollBar::valueChanged,
//            this,                  &DynamicView::notifyViewChanged);
//    connect(verticalScrollBar(), &QScrollBar::valueChanged,
//            this,                &DynamicView::notifyViewChanged);

    // if (QGLFormat::hasOpenGL())
    //    setViewport(new QGLWidget(QGLFormat(QGL::SampleBuffers)));

}

void DynamicView::wheelEvent(QWheelEvent* event)
{
    if (event->delta() > 0)
        changeZoom(+1);
    else
        changeZoom(-1);
}

void DynamicView::mouseDoubleClickEvent(QMouseEvent*)
{
    if (nextFit)
        fitZoom();
    else
        fullZoom();
}

void DynamicView::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton)
        mousePressed = true;
    QGraphicsView::mousePressEvent(event);
}

void DynamicView::mouseMoveEvent(QMouseEvent* event)
{
    QGraphicsView::mouseMoveEvent(event);
    if (mousePressed)
        notifyViewChanged();
}

void DynamicView::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton)
        mousePressed = false;
    QGraphicsView::mouseReleaseEvent(event);
}

void DynamicView::resizeEvent(QResizeEvent* event)
{
//    qDebug() << "Resize event";
    QGraphicsView::resizeEvent(event);
    if (matrix().m11() <= fitFactor)
        fitZoom();
    else
        notifyViewChanged();
}

// -------------------------------------------------------------------------

void DynamicView::fullZoom()
{
//    qDebug() << "Full zoom";
    changeZoom(0);
    nextFit = true;
}

void DynamicView::fitZoom()
{
//    qDebug() << "Fit zoom";
    fitInView(scene()->sceneRect(), Qt::KeepAspectRatio);
    fitFactor = matrix().m11();
    nextFit = false;
    notifyViewChanged();
}

void DynamicView::changeZoom(int direction)
{
//    qDebug() << "Change zoom";
    double factor = 1;
    if (direction != 0)
    {
        const double ZOOM_STEP = 0.25;
        double level = log2(matrix().m11());
        level += direction > 0 ? +ZOOM_STEP : -ZOOM_STEP;
        factor = pow(2, level);
        if (factor < fitFactor)
            factor = fitFactor;
    }
    setZoomFactor(factor);
    notifyViewChanged();
}

void DynamicView::setZoomFactor(double factor)
{
//    qDebug() << "Set zoom";
    QMatrix matrix;
    matrix.scale(factor, factor);
    setMatrix(matrix);
}

// -----------------------------------------------------------------------------

void DynamicView::changeView(QRect, double zoomFactor, int horizScroll, int vertScroll)
{
//    qDebug() << "Change view";
    double oldZoom = matrix().m11();
    int oldHoriz = horizontalScrollBar()->value();
    int oldVert = verticalScrollBar()->value();
    if (zoomFactor != oldZoom || horizScroll != oldHoriz || vertScroll != oldVert)
    {
        setZoomFactor(zoomFactor);
        horizontalScrollBar()->setValue(horizScroll);
        verticalScrollBar()->setValue(vertScroll);
        notifyViewChanged();
    }
}

void DynamicView::notifyViewChanged()
{
//    qDebug() << "Notify view";
    QRect sceneRect = getSceneRect();
    if (sceneRect.isEmpty())
    {
//        qDebug() << "NOT VALID";
        return;
    }
    int horizScroll = horizontalScrollBar()->value();
    int vertScroll = verticalScrollBar()->value();
    double zoomFactor = matrix().m11();
//    qDebug() << "COUNT =" << ++count << ", RECT =" << sceneRect << ", HORIZ =" << horizScroll << ", VERT =" << vertScroll << ", OK =" << !sceneRect.isEmpty();
    emit viewChanged(sceneRect, zoomFactor, horizScroll, vertScroll);
}

QRect DynamicView::getSceneRect()
{
//    qDebug() << "Get scene rect";
    QPoint topLeft = mapToScene(0, 0).toPoint();
    if (topLeft.x() < 0)
        topLeft.setX(0);
    if (topLeft.y() < 0)
        topLeft.setY(0);
    QSize viewSize = viewport()->size();
    QPoint bottomRight = mapToScene(viewSize.width(), viewSize.height()).toPoint();
    QRect imageSize = sceneRect().toRect();
    if (bottomRight.x() >= imageSize.width())
        bottomRight.setX(imageSize.width() - 1);
    if (bottomRight.y() >= imageSize.height())
        bottomRight.setY(imageSize.height() - 1);
    //qDebug() << sceneRect() << " ---" << QRect(topLeft, bottomRight);
    //qDebug() << "TL =" << topLeft << ", BR =" << bottomRight;
    //qDebug() << transform();
    QRect rect(topLeft, bottomRight);
    //if (rect.isEmpty())
    //    return imageSize;
        //qDebug() << rect;
    return rect;
}

/*
double DynamicView::getZoomFactor() const
{
    return matrix().m11();
}
*/

/*
double DynamicView::getZoomLevel() const
{
    return (100 * matrix().m11());
}
*/

/*
QRect DynamicView::getViewRect() const
{
    QPoint topLeft = mapToScene(0, 0).toPoint();
    if (topLeft.x() < 0)
        topLeft.setX(0);
    if (topLeft.y() < 0)
        topLeft.setY(0);
    QSize viewSize = viewport()->size();
    QPoint bottomRight = mapToScene(viewSize.width(), viewSize.height()).toPoint();
    QRect imageSize = sceneRect().toRect();
    if (bottomRight.x() >= imageSize.width())
        bottomRight.setX(imageSize.width() - 1);
    if (bottomRight.y() >= imageSize.height())
        bottomRight.setY(imageSize.height() - 1);
    //qDebug() << sceneRect() << " ---" << QRect(topLeft, bottomRight);
    //qDebug() << "TL =" << topLeft << ", BR =" << bottomRight;
    //qDebug() << transform();
    return QRect(topLeft, bottomRight);
}

void DynamicView::setViewRect(QRect rect)
{
    //QPolygon polygon = mapFromScene(rect);
    //QRect rect2(polygon.at(0), polygon.at(3));
    qDebug() << mapFromScene(rect);
}
*/
