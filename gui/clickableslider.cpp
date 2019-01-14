#include "clickableslider.h"

ClickableSlider::ClickableSlider(int min, int max, int ticks, int step1, int step2, QWidget* parent)
    : QSlider(Qt::Horizontal, parent)
{
    setRange(min, max);
    setTickPosition(QSlider::TicksBelow);
    setTickInterval((max - min + 1) / ticks);
    setSingleStep(step1);
    setPageStep(step2);
}

void ClickableSlider::mouseDoubleClickEvent(QMouseEvent *)
{
    emit doubleClicked();
}
