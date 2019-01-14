#include "echowidget.h"
#include "utility.h"

EchoWidget::EchoWidget(const cv::Mat &input, QWidget* parent)
    : ToolWidget(parent)
{
    radiusSpin = new QSpinBox();
    radiusSpin->setRange(1, 31);
    radiusSpin->setSingleStep(2);
    radiusSpin->setSuffix(tr(" px"));
    radiusSpin->setValue(5);
    connect(radiusSpin, static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged),
            this,       &EchoWidget::updateView);

    contrastSlider = new ClickableSlider(0, 255, 16, 1, 8);
    contrastSlider->setValue(215);
    connect(contrastSlider, &ClickableSlider::valueChanged,
            this,           &EchoWidget::updateView);
    contrastLabel = new QLabel(QString::number(contrastSlider->value()));

    QHBoxLayout* horizLayout = new QHBoxLayout();
    horizLayout->addWidget(new QLabel(tr("Radius:")));
    horizLayout->addWidget(radiusSpin);
    horizLayout->addWidget(new QLabel(tr("Contrast:")));
    horizLayout->addWidget(contrastSlider);
    horizLayout->addWidget(contrastLabel);
    horizLayout->addStretch();

    cv::Mat gray(input.rows, input.cols, CV_8UC3, cv::Scalar(127, 127, 127));
    echoViewer = new SingleViewer(input, gray);
    cv::split(input, channels);
    updateView();

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(horizLayout);
    vertLayout->addWidget(echoViewer);
    setLayout(vertLayout);
}

void EchoWidget::updateView()
{
    echoViewer->updateProcessed(computeEcho());
}

cv::Mat EchoWidget::computeEcho()
{
    QElapsedTimer timer;
    timer.start();

    int radius = radiusSpin->value();
    int contrast = contrastSlider->value();
    std::vector<cv::Mat> laplacians;
    cv::Mat lut = Utility::createLUT(0, contrast);
    for (uint i = 0; i < channels.size(); ++i)
    {
        cv::Mat deriv;
        cv::Laplacian(channels.at(i), deriv, CV_32F, radius);
        cv::normalize(abs(deriv), deriv, 0, 255, cv::NORM_MINMAX);
        deriv.convertTo(deriv, CV_8U);
        cv::LUT(deriv, lut, deriv);
        laplacians.push_back(deriv.clone());
    }
    cv::Mat edges;
    cv::merge(laplacians, edges);

    emit messageToShow(tr("[LUMINANCE::Echo] Processing time = %1 ms").arg(timer.elapsed()));
    return edges;
}
