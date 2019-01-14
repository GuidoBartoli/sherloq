#ifndef TOOLWIDGET_H
#define TOOLWIDGET_H

#include <QtWidgets>

class ToolWidget : public QWidget
{
    Q_OBJECT

public:
    explicit ToolWidget(QWidget* parent = Q_NULLPTR);

signals:
    void messageToShow(QString message);

};

#endif // TOOLWIDGET_H
