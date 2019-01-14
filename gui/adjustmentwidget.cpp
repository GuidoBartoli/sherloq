#include "adjustmentwidget.h"
#include "utility.h"

AdjustmentWidget::AdjustmentWidget(const cv::Mat &input, QWidget* parent)
    : ToolWidget(parent),
      input(input)
{
    cv::cvtColor(input, hsv, CV_BGR2HSV);
    brightnessSlider = new ClickableSlider(-255, +255, 32, 1, 8);
    brightnessLabel  = new QLabel(tr("0"));
    contrastSlider   = new ClickableSlider(-127, +127, 16, 1, 4);
    contrastLabel    = new QLabel(tr("0"));
    hueSlider        = new ClickableSlider(0, 180, 10, 1, 5);
    hueLabel         = new QLabel(tr("0°"));
    saturationSlider = new ClickableSlider(-255, +255, 32, 1, 8);
    saturationLabel  = new QLabel(tr("0"));
    invertCheck      = new QCheckBox(tr("Invert"));
    resetButton      = new QPushButton(tr("Reset"));
    resetSliders();

    connect(brightnessSlider, &QSlider::valueChanged,
            this,             &AdjustmentWidget::updateView);
    connect(brightnessSlider, &ClickableSlider::doubleClicked,
            this,             &AdjustmentWidget::resetBrightness);
    connect(contrastSlider,   &QSlider::valueChanged,
            this,             &AdjustmentWidget::updateView);
    connect(contrastSlider,   &ClickableSlider::doubleClicked,
            this,             &AdjustmentWidget::resetContrast);
    connect(hueSlider,        &QSlider::valueChanged,
            this,             &AdjustmentWidget::updateView);
    connect(hueSlider,        &ClickableSlider::doubleClicked,
            this,             &AdjustmentWidget::resetHue);
    connect(saturationSlider, &QSlider::valueChanged,
            this,             &AdjustmentWidget::updateView);
    connect(saturationSlider, &ClickableSlider::doubleClicked,
            this,             &AdjustmentWidget::resetSaturation);
    connect(invertCheck,      &QCheckBox::toggled,
            this,             &AdjustmentWidget::updateView);
    connect(resetButton,      &QPushButton::clicked,
            this,             &AdjustmentWidget::resetSliders);

    QGridLayout* gridLayout = new QGridLayout();
    gridLayout->addWidget(new QLabel(tr("Brightness:")), 0, 0);
    gridLayout->addWidget(brightnessSlider, 0, 1);
    gridLayout->addWidget(brightnessLabel, 0, 2);
    gridLayout->addWidget(new QLabel(tr("Contrast:")), 0, 3);
    gridLayout->addWidget(contrastSlider, 0, 4);
    gridLayout->addWidget(contrastLabel, 0, 5);
    gridLayout->addWidget(new QLabel(tr("Hue:")), 1, 0);
    gridLayout->addWidget(hueSlider, 1, 1);
    gridLayout->addWidget(hueLabel, 1, 2);
    gridLayout->addWidget(new QLabel(tr("Saturation:")), 1, 3);
    gridLayout->addWidget(saturationSlider, 1, 4);
    gridLayout->addWidget(saturationLabel, 1, 5);
    gridLayout->addWidget(invertCheck, 0, 6);
    gridLayout->addWidget(resetButton, 1, 6);

    adjustmentViewer = new SingleViewer(input, input);
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(gridLayout);
    vertLayout->addWidget(adjustmentViewer);
    setLayout(vertLayout);
}

void AdjustmentWidget::resetSliders()
{
    resetHue();
    resetBrightness();
    resetContrast();
    resetSaturation();
    invertCheck->setChecked(false);
}

void AdjustmentWidget::resetBrightness()
{
    brightnessSlider->setValue(0);
}

void AdjustmentWidget::resetContrast()
{
    contrastSlider->setValue(0);
}

void AdjustmentWidget::resetHue()
{
    hueSlider->setValue(0);
}

void AdjustmentWidget::resetSaturation()
{
    saturationSlider->setValue(0);
}

void AdjustmentWidget::updateView()
{
    QElapsedTimer timer;
    timer.start();

    int brightness = brightnessSlider->value();
    brightnessLabel->setText(tr("%1%2").arg(brightness > 0 ? "+" : "").arg(brightness));
    int contrast = contrastSlider->value();
    contrastLabel->setText(tr("%1%2").arg(contrast > 0 ? "+" : "").arg(contrast));
    int hue = hueSlider->value();
    hueLabel->setText(tr("%1°").arg(hue*2));
    int saturation = saturationSlider->value();
    saturationLabel->setText(tr("%1%2").arg(saturation > 0 ? "+" : "").arg(saturation));
    bool invert = invertCheck->isChecked();

    cv::Mat result;
    if (brightness != 0 || saturation != 0 || hue != 0)
    {
        std::vector<cv::Mat> chans;
        cv::split(hsv, chans);
        if (hue != 0)
        {
            chans[0].convertTo(chans[0], CV_32F);
            chans[0] += hue;
            cv::Mat mask = chans.at(0) < 0;
            cv::add(chans[0], +180, chans[0], mask);
            mask = chans.at(0) > 180;
            cv::add(chans[0], -180, chans[0], mask);
            chans[0].convertTo(chans[0], CV_8U);
        }
        if (saturation != 0)
            chans[1] += saturation;
        if (brightness != 0)
            chans[2] += brightness;
        cv::merge(chans, result);
        cv::cvtColor(result, result, CV_HSV2BGR);
    }
    else
        result = input.clone();
    if (contrast != 0)
    {
        cv::Mat lut = Utility::createLUT(contrast, contrast);
        cv::LUT(result, lut, result);
    }
    if (invert)
        cv::bitwise_not(result, result);

    emit messageToShow(tr("[INSPECTION::Adjustments] Processing time = %1 ms").arg(timer.elapsed()));
    adjustmentViewer->updateProcessed(result);
}
