#ifndef CLICKABLESLIDER_H
#define CLICKABLESLIDER_H

#include <QSlider>

class ClickableSlider : public QSlider
{
    Q_OBJECT

public:
    explicit ClickableSlider(int min = 0, int max = 100, int ticks = 10,
                             int step1 = 1, int step2 = 5, QWidget* parent = Q_NULLPTR);

protected:
    void mouseDoubleClickEvent(QMouseEvent*);

signals:
    void doubleClicked();
};

#endif // CLICKABLESLIDER_H
