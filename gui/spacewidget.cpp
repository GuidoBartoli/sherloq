#include "spacewidget.h"

SpaceWidget::SpaceWidget(const cv::Mat &input, QWidget* parent)
    : ToolWidget(parent)
{
    cv::Mat input0, conv;
    input.convertTo(input0, CV_32FC3, 1 / 255.0);

    cv::split(input, rgb);
    std::swap(rgb[0], rgb[2]);

    cv::cvtColor(input0, conv, CV_BGR2HSV);
    cv::split(conv, hsv);
    hsv[0].convertTo(hsv[0], CV_8U, 255.0 / 360.0);
    hsv[1].convertTo(hsv[1], CV_8U, 255);
    hsv[2].convertTo(hsv[2], CV_8U, 255);

    cv::cvtColor(input0, conv, CV_BGR2HLS);
    cv::split(conv, hls);
    hls[0].convertTo(hls[0], CV_8U, 255.0 / 360.0);
    hls[1].convertTo(hls[1], CV_8U, 255);
    hls[2].convertTo(hls[2], CV_8U, 255);

    for (int i = 0; i < 4; ++i)
        cmyk.push_back(cv::Mat(input.size(), CV_8U));

    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            float r = (int)rgb.at(0).at<uchar>(i, j) / 255.0;
            float g = (int)rgb.at(1).at<uchar>(i, j) / 255.0;
            float b = (int)rgb.at(2).at<uchar>(i, j) / 255.0;
            float k = std::min(std::min(1 - r, 1 - g), 1 - b);
            cmyk[0].at<uchar>(i, j) = (1 - r - k) / (1 - k) * 255.0;
            cmyk[1].at<uchar>(i, j) = (1 - g - k) / (1 - k) * 255.0;
            cmyk[2].at<uchar>(i, j) = (1 - b - k) / (1 - k) * 255.0;
            cmyk[3].at<uchar>(i, j) = k * 255.0;
        }
    }

    cv::cvtColor(input, conv, CV_BGR2YCrCb);
    cv::split(conv, ycrcb);

    cv::cvtColor(input, conv, CV_BGR2XYZ);
    cv::split(conv, xyz);

    cv::cvtColor(input, conv, CV_BGR2Lab);
    cv::split(conv, lab);

    cv::cvtColor(input, conv, CV_BGR2Luv);
    cv::split(conv, luv);

    rgbRadio = new QRadioButton(tr("Red-Green-Blue"));
    rgbCombo = new QComboBox();
    rgbCombo->addItem(tr("R"));
    rgbCombo->addItem(tr("G"));
    rgbCombo->addItem(tr("B"));
    rgbRadio->setChecked(true);
    lastRadio = rgbRadio;
    connect(rgbRadio, &QRadioButton::toggled,
            this,     &SpaceWidget::updateView);
    connect(rgbCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this,     &SpaceWidget::updateView);

    cmykRadio = new QRadioButton(tr("CMYK Inks"));
    cmykCombo = new QComboBox();
    cmykCombo->addItem(tr("C"));
    cmykCombo->addItem(tr("M"));
    cmykCombo->addItem(tr("Y"));
    cmykCombo->addItem(tr("K"));
    connect(cmykRadio, &QRadioButton::toggled,
            this,      &SpaceWidget::updateView);
    connect(cmykCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this,      &SpaceWidget::updateView);

    grayRadio = new QRadioButton(tr("Grayscale"));
    grayCombo = new QComboBox();
    grayCombo->addItem(tr("Light"));
    grayCombo->addItem(tr("Luma"));
    grayCombo->addItem(tr("Mean"));

    hsvRadio = new QRadioButton(tr("Hue-Sat-Val"));
    hsvCombo = new QComboBox();
    hsvCombo->addItem(tr("H"));
    hsvCombo->addItem(tr("S"));
    hsvCombo->addItem(tr("V"));
    connect(hsvRadio, &QRadioButton::toggled,
            this,     &SpaceWidget::updateView);
    connect(hsvCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this,     &SpaceWidget::updateView);

    hlsRadio = new QRadioButton(tr("Hue-Lum-Sat"));
    hlsCombo = new QComboBox();
    hlsCombo->addItem(tr("H"));
    hlsCombo->addItem(tr("L"));
    hlsCombo->addItem(tr("S"));
    connect(hlsRadio, &QRadioButton::toggled,
            this,     &SpaceWidget::updateView);
    connect(hlsCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this,     &SpaceWidget::updateView);

    ycrcbRadio = new QRadioButton(tr("YCrCb JPEG"));
    ycrcbCombo = new QComboBox();
    ycrcbCombo->addItem(tr("Y"));
    ycrcbCombo->addItem(tr("Cr"));
    ycrcbCombo->addItem(tr("Cb"));
    connect(ycrcbRadio, &QRadioButton::toggled,
            this,       &SpaceWidget::updateView);
    connect(ycrcbCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this,       &SpaceWidget::updateView);

    xyzRadio = new QRadioButton(tr("CIE XYZ.Rec709"));
    xyzCombo = new QComboBox();
    xyzCombo->addItem(tr("X"));
    xyzCombo->addItem(tr("Y"));
    xyzCombo->addItem(tr("Z"));
    connect(xyzRadio, &QRadioButton::toggled,
            this,     &SpaceWidget::updateView);
    connect(xyzCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this,     &SpaceWidget::updateView);

    labRadio = new QRadioButton(tr("CIE L*a*b*"));
    labCombo = new QComboBox();
    labCombo->addItem(tr("L"));
    labCombo->addItem(tr("a"));
    labCombo->addItem(tr("b"));
    connect(labRadio, &QRadioButton::toggled,
            this,     &SpaceWidget::updateView);
    connect(labCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this,     &SpaceWidget::updateView);

    luvRadio = new QRadioButton(tr("CIE L*u*v*"));
    luvCombo = new QComboBox();
    luvCombo->addItem(tr("L"));
    luvCombo->addItem(tr("u"));
    luvCombo->addItem(tr("v"));
    connect(luvRadio, &QRadioButton::toggled,
            this,     &SpaceWidget::updateView);
    connect(luvCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this,     &SpaceWidget::updateView);

    QGridLayout* gridLayout = new QGridLayout();
    gridLayout->addWidget(rgbRadio, 0, 0);
    gridLayout->addWidget(rgbCombo, 0, 1);
    gridLayout->addWidget(cmykRadio, 0, 2);
    gridLayout->addWidget(cmykCombo, 0, 3);
    gridLayout->addWidget(grayRadio, 0, 4);
    gridLayout->addWidget(grayCombo, 0, 5);
    gridLayout->addWidget(hsvRadio, 1, 0);
    gridLayout->addWidget(hsvCombo, 1, 1);
    gridLayout->addWidget(hlsRadio, 1, 2);
    gridLayout->addWidget(hlsCombo, 1, 3);
    gridLayout->addWidget(ycrcbRadio, 1, 4);
    gridLayout->addWidget(ycrcbCombo, 1, 5);
    gridLayout->addWidget(xyzRadio, 2, 0);
    gridLayout->addWidget(xyzCombo, 2, 1);
    gridLayout->addWidget(labRadio, 2, 2);
    gridLayout->addWidget(labCombo, 2, 3);
    gridLayout->addWidget(luvRadio, 2, 4);
    gridLayout->addWidget(luvCombo, 2, 5);

    colorViewer = new SingleViewer(input, rgb.at(0));

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(gridLayout);
    vertLayout->addWidget(colorViewer);
    setLayout(vertLayout);
}

void SpaceWidget::updateView()
{
    if (rgbRadio->isChecked())
    {
        colorViewer->updateProcessed(rgb.at(rgbCombo->currentIndex()));
        lastRadio = rgbRadio;
    }
    else if (cmykRadio->isChecked())
    {
        colorViewer->updateProcessed(cmyk.at(cmykCombo->currentIndex()));
        lastRadio = cmykRadio;
    }
    else if (hsvRadio->isChecked())
    {
        colorViewer->updateProcessed(hsv.at(hsvCombo->currentIndex()));
        lastRadio = hsvRadio;
    }
    else if (hlsRadio->isChecked())
    {
        colorViewer->updateProcessed(hls.at(hlsCombo->currentIndex()));
        lastRadio = hlsRadio;
    }
    else if (ycrcbRadio->isChecked())
    {
        colorViewer->updateProcessed(ycrcb.at(ycrcbCombo->currentIndex()));
        lastRadio = ycrcbRadio;
    }
    else if (luvRadio->isChecked())
    {
        colorViewer->updateProcessed(luv.at(luvCombo->currentIndex()));
        lastRadio = luvRadio;
    }
    else if (xyzRadio->isChecked())
    {
        colorViewer->updateProcessed(xyz.at(xyzCombo->currentIndex()));
        lastRadio = xyzRadio;
    }
    else if (labRadio->isChecked())
    {
        colorViewer->updateProcessed(lab.at(labCombo->currentIndex()));
        lastRadio = labRadio;
    }
    else
        lastRadio->setChecked(true);
}
