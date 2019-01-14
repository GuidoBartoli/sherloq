#include "elawidget.h"
#include "utility.h"

ElaWidget::ElaWidget(const cv::Mat &input, QWidget* parent)
    : ToolWidget(parent),
      input(input)
{
    qualitySpin = new QSpinBox();
    qualitySpin->setRange(0, 100);
    scaleSpin = new QSpinBox();
    scaleSpin->setRange(0, 100);
    resetParams();
    QPushButton* resetButton = new QPushButton(tr("Default"));
    connect(resetButton, &QPushButton::clicked,
            this,        &ElaWidget::resetParams);

    QHBoxLayout* paramsLayout = new QHBoxLayout();
    paramsLayout->addWidget(new QLabel(tr("JPEG Quality:")));
    paramsLayout->addWidget(qualitySpin);
    paramsLayout->addWidget(new QLabel(tr("Error Scale:")));
    paramsLayout->addWidget(scaleSpin);
    paramsLayout->addWidget(resetButton);
    paramsLayout->addStretch();

    elaViewer = new SingleViewer(input, computeEla());

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(paramsLayout);
    vertLayout->addWidget(elaViewer);
    setLayout(vertLayout);

    connect(qualitySpin, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged),
            this,        &ElaWidget::updateView);
    connect(scaleSpin,   static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged),
            this,        &ElaWidget::updateView);
}

void ElaWidget::updateView(int)
{
    elaViewer->updateProcessed(computeEla());
}

cv::Mat ElaWidget::computeEla()
{
    QElapsedTimer timer;
    timer.start();

    int quality = qualitySpin->value();
    int scale = scaleSpin->value();
    cv::Mat compressed = Utility::jpegCompress(input, quality);
    cv::Mat ela = cv::abs(compressed - input) * scale;

    emit messageToShow(tr("[JPEG::ELA] Processing time = %1 ms").arg(timer.elapsed()));
    return ela;
}

void ElaWidget::resetParams()
{
    qualitySpin->setValue(DEFAULT_QUALITY);
    scaleSpin->setValue(DEFAULT_SCALE);
}
