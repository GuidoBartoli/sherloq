#include "sweepwidget.h"
#include "utility.h"

SweepWidget::SweepWidget(cv::Mat &input, QWidget *parent)
    : ToolWidget(parent),
      input(input)
{
    sweepSlider = new QSlider(Qt::Horizontal);
    sweepSlider->setRange(0, 100);
    sweepSlider->setSingleStep(1);
    sweepSlider->setPageStep(4);
    sweepSlider->setTickPosition(QSlider::TicksBelow);
    sweepSlider->setTickInterval(10);
    sweepLabel = new QLabel();
    widthSpin = new QSpinBox();
    widthSpin->setRange(8, 128);
    resetParams();
    QPushButton* resetButton = new QPushButton(tr("Default"));
    connect(resetButton, &QPushButton::clicked,
            this,        &SweepWidget::resetParams);

    QHBoxLayout* horizLayout = new QHBoxLayout();
    horizLayout->addWidget(new QLabel(tr("Sweep:")));
    horizLayout->addWidget(sweepSlider);
    horizLayout->addWidget(sweepLabel);
    horizLayout->addWidget(new QLabel(tr("Width:")));
    horizLayout->addWidget(widthSpin);
    horizLayout->addWidget(resetButton);
    horizLayout->addStretch();

    sweepViewer = new SingleViewer(input, computeSweep());

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(horizLayout);
    vertLayout->addWidget(sweepViewer);
    setLayout(vertLayout);

    connect(sweepSlider, &QSlider::valueChanged,
            this,        &SweepWidget::updateView);
    connect(widthSpin, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged),
            this,      &SweepWidget::updateView);
}

cv::Mat SweepWidget::computeSweep()
{
    QElapsedTimer timer;
    timer.start();
    uchar mid = (uchar)(sweepSlider->value() / 100.0 * 255);
    int width = widthSpin->value();
    sweepLabel->setText(tr("%1%").arg(sweepSlider->value()));
    uchar low = cv::saturate_cast<uchar>(mid - width);
    uchar high = cv::saturate_cast<uchar>(mid + width);
    cv::Mat lut = Utility::createLUT(low, 255 - high);
    cv::Mat contrast;
    cv::LUT(input, lut, contrast);
    emit messageToShow(tr("[INSPECTION::Sweep] Processing time = %1 ms").arg(timer.elapsed()));
    return contrast;
}

void SweepWidget::updateView(int)
{
    sweepViewer->updateProcessed(computeSweep());
}

void SweepWidget::resetParams()
{
    sweepSlider->setValue(50);
    widthSpin->setValue(32);
}
