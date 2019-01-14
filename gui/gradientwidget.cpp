#include "gradientwidget.h"
#include "utility.h"

GradientWidget::GradientWidget(const cv::Mat &input, QWidget* parent)
    : ToolWidget(parent),
      input(input)
{
    intensitySpin = new QSpinBox();
    intensitySpin->setRange(-1, 16);
    intensitySpin->setSuffix(tr(" levels"));
    intensitySpin->setSpecialValueText(tr("Equalized"));
    invertCheck = new QCheckBox(tr("Invert"));
    modulusCheck = new QCheckBox(tr("Modulus"));
    //equalizeCheck = new QCheckBox(tr("Equalize"));
    QPushButton* resetButton = new QPushButton(tr("Default"));
    resetParams();

    QHBoxLayout* horizLayout = new QHBoxLayout();
    horizLayout->addWidget(new QLabel(tr("Intensity:")));
    horizLayout->addWidget(intensitySpin);
    horizLayout->addWidget(invertCheck);
    horizLayout->addWidget(modulusCheck);
    //horizLayout->addWidget(equalizeCheck);
    horizLayout->addStretch();

    gradientViewer = new SingleViewer(input, computeGradient());
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(horizLayout);
    vertLayout->addWidget(gradientViewer);
    setLayout(vertLayout);

    connect(intensitySpin,  static_cast<void(QSpinBox::*)(int)>(&QSpinBox::valueChanged),
            this,           &GradientWidget::updateView);
    connect(invertCheck,    &QCheckBox::toggled,
            this,           &GradientWidget::updateView);
    connect(modulusCheck,   &QCheckBox::toggled,
            this,           &GradientWidget::updateView);
    //connect(equalizeCheck,  &QCheckBox::toggled,
    //        this,           &GradientWidget::updateView);
    connect(resetButton,    &QPushButton::clicked,
            this,           &GradientWidget::resetParams);
}

void GradientWidget::updateView(int)
{
    gradientViewer->updateProcessed(computeGradient());
}

void GradientWidget::resetParams()
{
    intensitySpin->setValue(-1);
    invertCheck->setChecked(false);
    modulusCheck->setChecked(false);
    //equalizeCheck->setChecked(true);
}

cv::Mat GradientWidget::computeGradient()
{
    QElapsedTimer timer;
    timer.start();

    int intensity = intensitySpin->value();
    bool invert = invertCheck->isChecked();
    bool modulus = modulusCheck->isChecked();
    intensitySpin->setEnabled(!modulus);
    invertCheck->setEnabled(!modulus);
    //bool equalize = equalizeCheck->isChecked();

    if (scharrX.empty() || scharrY.empty() || scharrZ.empty())
    {
        cv::Mat gray;
        cv::cvtColor(input, gray, CV_BGR2GRAY);
        cv::Scharr(gray, scharrX, CV_32F, 1, 0);
        cv::Scharr(gray, scharrY, CV_32F, 0, 1);
        scharrZ = abs(scharrX + scharrY);
        cv::normalize(scharrZ, scharrZ, 0, 255, cv::NORM_MINMAX);
        scharrZ.convertTo(scharrZ, CV_8U);
        cv::Mat lutZ = Utility::createLUT();
        cv::LUT(scharrZ, lutZ, scharrZ);
    }

    cv::Mat dX = scharrX * (!invert ? -1 : +1);
    cv::Mat dY = scharrY * (!invert ? -1 : +1);
    cv::Mat dZ = scharrZ;

    if (!modulus)
    {
        double min, max;
        cv::minMaxLoc(dX, &min, &max);
        double limit = std::max(fabs(min), fabs(max));
        dX = (dX / limit + 1) * 127;
        cv::minMaxLoc(dY, &min, &max);
        limit = std::max(fabs(min), fabs(max));
        dY = (dY / limit + 1) * 127;
    }
    else
    {
        dX = abs(dX.clone());
        dY = abs(dY.clone());
        cv::normalize(dX, dX, 0, 255, cv::NORM_MINMAX);
        cv::normalize(dY, dY, 0, 255, cv::NORM_MINMAX);
    }
    dX.convertTo(dX, CV_8U);
    dY.convertTo(dY, CV_8U);

    if (intensity >= 0 && !modulus)
    {
        int maxIntensity = intensitySpin->maximum();
        uchar low = 127 - (maxIntensity - intensity);
        uchar high = 127 + (maxIntensity - intensity);
        cv::Mat lutXY = Utility::createLUT(low, 255 - high);
        cv::LUT(dX, lutXY, dX);
        cv::LUT(dY, lutXY, dY);
    }
    else
    {
        cv::equalizeHist(dX, dX);
        cv::equalizeHist(dY, dY);
    }

    cv::Mat output;
    std::vector<cv::Mat> gradient;
    gradient.push_back(dZ);
    //gradient.push_back(!invert ? dX : dY);
    //gradient.push_back(!invert ? dY : dX);
    gradient.push_back(dX);
    gradient.push_back(dY);
    cv::merge(gradient, output);

    emit messageToShow(tr("[LUMINANCE::Gradient] Processing time = %1 ms").arg(timer.elapsed()));
    return output;

//    cv::imwrite("/home/bartoli/Downloads/dx.png", dx);
//    cv::imwrite("/home/bartoli/Downloads/dy.png", dy);





    /*
    cv::Mat gray;
    cv::cvtColor(input, gray, CV_BGR2GRAY);
    cv::Mat kernel = (cv::Mat_<float>(3, 3) <<
                       3,  10,  3,
                       0,   0,  0,
                      -3, -10, -3);
    cv::Mat dx, dy;
    cv::filter2D(gray, dx, CV_32F, kernel);
    cv::filter2D(gray, dy, CV_32F, kernel.t());
    double min, max;
    cv::minMaxLoc(dx, &min, &max);
    qDebug() << min << max;
    cv::normalize(dx, dx, 0, 255, cv::NORM_MINMAX);
    cv::normalize(dy, dy, 0, 255, cv::NORM_MINMAX);
    dx.convertTo(dx, CV_8U);
    dy.convertTo(dy, CV_8U);
    cv::imwrite("/home/bartoli/Downloads/dx.png", dx);
    cv::imwrite("/home/bartoli/Downloads/dy.png", dy);
    return dx;
    */

    /*
    cv::Mat gray;
    cv::cvtColor(input, gray, CV_BGR2GRAY);
    cv::Mat dx, dy;
    cv::Scharr(gray, dx, CV_32F, 1, 0);
    cv::Scharr(gray, dy, CV_32F, 0, 1);
    cv::normalize(dx, dx, 0, 255, cv::NORM_MINMAX);
    cv::normalize(dy, dy, 0, 255, cv::NORM_MINMAX);
    dx.convertTo(dx, CV_8U);
    dy.convertTo(dy, CV_8U);
    cv::imwrite("/home/bartoli/Downloads/dx.png", dx);
    cv::imwrite("/home/bartoli/Downloads/dy.png", dy);
    return dx;
    */

    /*
    int value = intensitySpin->value();
    int ksize;
    if (value == 0)
        ksize = CV_SCHARR;
    else if (value == 1)
        ksize = 1;
    else if (value == 2)
        ksize = 3;
    else if (value == 3)
        ksize = 5;
    else if (value == 4)
        ksize = 7;
    cv::Mat kernel1, kernel2;
    cv::getDerivKernels(kernel1, kernel2, 1, 1, ksize);

    cv::Mat gray;
    cv::cvtColor(input, gray, CV_BGR2GRAY);
    cv::Mat dx, dy;
    cv::filter2D(gray, dx, CV_32F, kernel1);
    cv::filter2D(gray, dy, CV_32F, kernel2);
    cv::normalize(dx, dx, 0, 255, cv::NORM_MINMAX);
    cv::normalize(dy, dy, 0, 255, cv::NORM_MINMAX);
    dx.convertTo(dx, CV_8U);
    dy.convertTo(dy, CV_8U);
    cv::imwrite("/home/bartoli/Downloads/dx.png", dx);
    cv::imwrite("/home/bartoli/Downloads/dy.png", dy);

    return dx;
    */

    /*
    const float COS2 = 0.707106781187;
    cv::Mat gray;
    cv::cvtColor(input, gray, CV_BGR2GRAY);
    gray.convertTo(gray, CV_32F, 1 / 255.0);

    cv::Mat c1(gray.rows, gray.cols, CV_32F);
    cv::Mat c2(gray.rows, gray.cols, CV_32F);
    cv::Mat c3(gray.rows, gray.cols, CV_32F);
    for (int i = 1; i < gray.rows - 1; ++i)
    {
        for (int j = 1; j < gray.cols - 1; ++j)
        {
            float x = 0, y = 0, z = 0;
            float tl = gray.at<float>(i - 1, j - 1);
            x -= tl * COS2;
            y -= tl * COS2;
            float t  = gray.at<float>(i - 1, j);
            y -= t;
            float tr = gray.at<float>(i - 1, j + 1);
            x += tr * COS2;
            y -= tr * COS2;
            float r  = gray.at<float>(i, j + 1);
            x += r;
            float br = gray.at<float>(i + 1, j + 1);
            x += br * COS2;
            y += br * COS2;
            float b  = gray.at<float>(i + 1, j);
            y += b;
            float bl = gray.at<float>(i + 1, j - 1);
            x -= bl * COS2;
            y += bl * COS2;
            float l  = gray.at<float>(i, j - 1);
            x -= l;
            z = sqrt(x*x + y*y);
            c1.at<float>(i, j) = -x;
            c2.at<float>(i, j) = -y;
            c3.at<float>(i, j) = z;
        }
    }
    double min, max;
    cv::minMaxLoc(c1, &min, &max);
    qDebug() << min << max;
    cv::normalize(c1, c1, 0, 255, cv::NORM_MINMAX);
    cv::normalize(c2, c2, 0, 255, cv::NORM_MINMAX);
    cv::normalize(c3, c3, 0, 255, cv::NORM_MINMAX);
    c1.convertTo(c1, CV_8U);
    c2.convertTo(c2, CV_8U);
    c3.convertTo(c3, CV_8U);
    //cv::equalizeHist(c1, c1);
    //cv::equalizeHist(c2, c2);
    c3 *= 2;
    //cv::equalizeHist(c3, c3);
    cv::imwrite("/tmp/c1.png", c1);
    cv::imwrite("/tmp/c2.png", c2);
    cv::imwrite("/tmp/c3.png", c3);

    std::vector<cv::Mat> gradient;
    gradient.push_back(c3);
    gradient.push_back(c1);
    gradient.push_back(c2);
    cv::Mat output;
    cv::merge(gradient, output);
    cv::imwrite("/home/bartoli/Downloads/gradient0.png", output);
    return output;
    */

    /*
    double intensity = 1 - (intensitySpin->value() / 100.0);
    bool invert = invertCheck->isChecked();
    bool normalize = normalizeCheck->isChecked();
    bool modulus = modulusCheck->isChecked();
    bool equalize = equalizeCheck->isChecked();

    cv::Mat gray;
    cv::cvtColor(input, gray, CV_BGR2GRAY);
    cv::Mat kernel1, kernel2;
    if (!invert)
    {
        kernel1 = (cv::Mat_<float>(5, 1) << 0, +1, 0, -1, 0);
        kernel2 = kernel1.t();
    }
    else
    {
        kernel2 = (cv::Mat_<float>(5, 1) << 0, -1, 0, +1, 0);
        kernel1 = kernel2.t();
    }
    uchar low = 127 - 127 * intensity;
    uchar high = 127 + 127 * intensity;
    cv::Mat lutXY = Utility::buildContrastLUT(low, high);

    cv::Mat red;
    cv::filter2D(gray, red, CV_32F, kernel1);
    if (modulus)
        red = cv::abs(red.clone());
    else
        red = (red + 255) / 2;
    if (normalize)
        cv::normalize(red, red, 0, 255, cv::NORM_MINMAX);
    red.convertTo(red, CV_8U);
    cv::LUT(red, lutXY, red);
    cv::imwrite("/tmp/red.png", red);
    if (equalize)
        cv::equalizeHist(red, red);

    cv::Mat green;
    cv::filter2D(gray, green, CV_32F, kernel2);
    if (modulus)
        green = cv::abs(green.clone());
    else
        green = (green + 255) / 2;
    if (normalize)
        cv::normalize(green, green, 0, 255, cv::NORM_MINMAX);
    green.convertTo(green, CV_8U);
    cv::LUT(green, lutXY, green);
    cv::imwrite("/tmp/green.png", green);
    if (equalize)
        cv::equalizeHist(green, green);

    cv::Mat blue1, blue2;
    cv::Sobel(gray, blue1, CV_32F, 1, 0, CV_SCHARR);
    cv::Sobel(gray, blue2, CV_32F, 0, 1, CV_SCHARR);
    cv::Mat blue = abs(blue1 + blue2);
    cv::normalize(blue, blue, 0, 255, cv::NORM_MINMAX);
    blue.convertTo(blue, CV_8U);
    cv::Mat lut2 = Utility::buildContrastLUT();
    cv::LUT(blue, lut2, blue);

    cv::Mat output;
    std::vector<cv::Mat> gradient;
    gradient.push_back(blue);
    gradient.push_back(green);
    gradient.push_back(red);
    cv::merge(gradient, output);

    //output *= intensity;

    return output;
    */

}
