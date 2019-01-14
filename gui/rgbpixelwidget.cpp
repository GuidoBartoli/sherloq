#include "rgbpixelwidget.h"

RgbPixelWidget::RgbPixelWidget(const cv::Mat &input, QWidget *parent)
    : ToolWidget(parent)
{
    minImage = cv::Mat::zeros(input.rows, input.cols, CV_8UC3);
    maxImage = cv::Mat::zeros(input.rows, input.cols, CV_8UC3);
    avgImage = cv::Mat::zeros(input.rows, input.cols, CV_8UC3);
    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            cv::Vec3b bgr = input.at<cv::Vec3b>(i, j);
            uchar b = bgr[0];
            uchar g = bgr[1];
            uchar r = bgr[2];

            if (b < g && b < r)
                minImage.at<cv::Vec3b>(i, j) = cv::Vec3b(b, 0, 0);
            else if (g < r && g < b)
                minImage.at<cv::Vec3b>(i, j) = cv::Vec3b(0, g, 0);
            else if (r < b && r < g)
                minImage.at<cv::Vec3b>(i, j) = cv::Vec3b(0, 0, r);

            if (b > g && b > r)
                maxImage.at<cv::Vec3b>(i, j) = cv::Vec3b(b, 0, 0);
            else if (g > r && g > b)
                maxImage.at<cv::Vec3b>(i, j) = cv::Vec3b(0, g, 0);
            else if (r > b && r > g)
                maxImage.at<cv::Vec3b>(i, j) = cv::Vec3b(0, 0, r);

            if ((b > r && b < g) || (b > g && b < r))
                avgImage.at<cv::Vec3b>(i, j) = cv::Vec3b(b, 0, 0);
            else if ((g > r && g < b) || (g > b && g < r))
                avgImage.at<cv::Vec3b>(i, j) = cv::Vec3b(0, g, 0);
            else if ((r > b && r < g) || (r > g && r < b))
                avgImage.at<cv::Vec3b>(i, j) = cv::Vec3b(0, 0, r);
        }
    }

    minButton = new QRadioButton(tr("Minimum"));
    minButton->setToolTip(tr("Pixels are colored as the RGB channel with minimum value"));
    maxButton = new QRadioButton(tr("Maximum"));
    maxButton->setToolTip(tr("Pixels are colored as the RGB channel with maximum value"));
    avgButton = new QRadioButton(tr("Average"));
    avgButton->setToolTip(tr("Pixels are colored as the RGB channel with average value"));
    minButton->setChecked(true);
    QHBoxLayout* horizLayout = new QHBoxLayout();
    horizLayout->addWidget(new QLabel(tr("Statistics:")));
    horizLayout->addWidget(minButton);
    horizLayout->addWidget(maxButton);
    horizLayout->addWidget(avgButton);
    horizLayout->addStretch();

    rgbPixelViewer = new SingleViewer(input, minImage);

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(horizLayout);
    vertLayout->addWidget(rgbPixelViewer);
    setLayout(vertLayout);

    connect(minButton, &QRadioButton::toggled,
            this,      &RgbPixelWidget::updateView);
    connect(maxButton, &QRadioButton::toggled,
            this,      &RgbPixelWidget::updateView);
    connect(avgButton, &QRadioButton::toggled,
            this,      &RgbPixelWidget::updateView);
}

void RgbPixelWidget::updateView(bool)
{
    if (minButton->isChecked())
    {
        rgbPixelViewer->updateProcessed(minImage);
        lastRadio = minButton;
    }
    else if (maxButton->isChecked())
    {
        rgbPixelViewer->updateProcessed(maxImage);
        lastRadio = maxButton;
    }
    else if (avgButton->isChecked())
    {
        rgbPixelViewer->updateProcessed(avgImage);
        lastRadio = avgButton;
    }
    else
        lastRadio->setChecked(true);
}
