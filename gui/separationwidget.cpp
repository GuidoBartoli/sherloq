#include "separationwidget.h"
#include "utility.h"

SeparationWidget::SeparationWidget(const cv::Mat &input, QWidget *parent)
    : ToolWidget(parent)
{
    cv::cvtColor(input, gray, CV_BGR2GRAY);
    cv::Mat denoised;
    cv::medianBlur(gray, denoised, 3);
    gray.convertTo(gray, CV_32F);
    denoised.convertTo(denoised, CV_32F);
    noise = cv::abs(gray - denoised);

    intensitySlider = new ClickableSlider(0, 255, 16, 1, 8);
    intensityLabel = new QLabel(tr("127"));
    accurateCheck = new QCheckBox(tr("Enhanced"));
    QPushButton* resetButton = new QPushButton(tr("Default"));
    resetParams();

    connect(intensitySlider, &QSlider::valueChanged,
            this,            &SeparationWidget::updateView);
    connect(intensitySlider, &ClickableSlider::doubleClicked,
            this,            &SeparationWidget::resetIntensity);
    connect(accurateCheck, &QCheckBox::toggled,
            this,          &SeparationWidget::updateView);
    connect(resetButton, &QPushButton::clicked,
            this,        &SeparationWidget::resetParams);

    QHBoxLayout* horizLayout = new QHBoxLayout();
    horizLayout->addWidget(new QLabel(tr("Intensity:")));
    horizLayout->addWidget(intensitySlider);
    horizLayout->addWidget(intensityLabel);
    horizLayout->addWidget(accurateCheck);
    horizLayout->addStretch();

    noiseViewer = new SingleViewer(input, computeNoise());

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(horizLayout);
    vertLayout->addWidget(noiseViewer);
    setLayout(vertLayout);
}

cv::Mat SeparationWidget::computeNoise()
{
    QElapsedTimer timer;
    timer.start();

    bool accurate = accurateCheck->isChecked();
    int intensity = accurate ? NOISE_CONT : intensitySlider->value();
    intensitySlider->setEnabled(!accurate);
    if (!accurate)
        intensityLabel->setNum(intensity);
    intensityLabel->setEnabled(!accurate);

    cv::Mat noise2 = noise.clone();
    if (accurate)
    {
        if (noise0.empty())
        {
            QProgressDialog progress(tr("Separating noise..."), QString(), 0, gray.rows, this);
            progress.setWindowModality(Qt::WindowModal);
            for (int i = BLOCK_PAD; i < gray.rows - BLOCK_PAD; ++i)
            {
                for (int j = BLOCK_PAD; j < gray.cols - BLOCK_PAD; ++j)
                {
                    cv::Mat block(gray, cv::Rect(j - BLOCK_PAD, i - BLOCK_PAD,
                                                 BLOCK_SIZE, BLOCK_SIZE));
                    cv::Scalar meanScalar, devScalar;
                    cv::meanStdDev(block, meanScalar, devScalar);
                    float stddev = devScalar[0] / MAX_STDDEV;
                    if (stddev != 0)
                        noise2.at<float>(i, j) /= stddev;
                    else
                        noise2.at<float>(i, j) = 0;
                }
                if (i % 100 == 0)
                    progress.setValue(i);
            }
            cv::GaussianBlur(noise2, noise2, cv::Size(BLOCK_SIZE, BLOCK_SIZE), 0);
            noise0 = noise2;
        }
        else
            noise2 = noise0.clone();
    }

    cv::normalize(noise2, noise2, 0, 255, cv::NORM_MINMAX);
    noise2.convertTo(noise2, CV_8U);
    cv::Mat lut = Utility::createLUT(0, intensity);
    cv::LUT(noise2, lut, noise2);

    emit messageToShow(tr("[NOISE::Separation] Processing time = %1 ms").arg(timer.elapsed()));
    return noise2;
}

void SeparationWidget::resetParams()
{
    accurateCheck->setChecked(false);
    resetIntensity();
}

void SeparationWidget::resetIntensity()
{
    intensitySlider->setValue(192);
}

void SeparationWidget::updateView()
{
    noiseViewer->updateProcessed(computeNoise());
}
