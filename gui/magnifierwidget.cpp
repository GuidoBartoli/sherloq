#include "magnifierwidget.h"
#include "utility.h"

MagnifierWidget::MagnifierWidget(const cv::Mat &input, QWidget *parent)
    : ToolWidget(parent),
      input(input)
{
    equalizeRadio = new QRadioButton(tr("Histogram Equalization"));
    contrast1Radio = new QRadioButton(tr("Auto Contrast (luminosity)"));
    contrast2Radio = new QRadioButton(tr("Auto Contrast (per channel)"));
    retinexRadio = new QRadioButton(tr("PDE-Retinex Algorithm"));
    equalizeRadio->setChecked(true);
    lastRadio = equalizeRadio;
    retinexRadio->setEnabled(false);

    QHBoxLayout* horizLayout = new QHBoxLayout();
    horizLayout->addWidget(equalizeRadio);
    horizLayout->addWidget(contrast1Radio);
    horizLayout->addWidget(contrast2Radio);
    horizLayout->addWidget(retinexRadio);
    horizLayout->addStretch();

    magnifierViewer = new SingleViewer(input, input);
    connect(magnifierViewer, &SingleViewer::viewChanged,
            this,            &MagnifierWidget::updateView);

    connect(equalizeRadio, &QRadioButton::toggled,
            this,          &MagnifierWidget::updateRect);
    connect(contrast1Radio, &QRadioButton::toggled,
            this,           &MagnifierWidget::updateRect);
    connect(contrast2Radio, &QRadioButton::toggled,
            this,           &MagnifierWidget::updateRect);
    connect(retinexRadio, &QRadioButton::toggled,
            this,         &MagnifierWidget::updateRect);

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(horizLayout);
    vertLayout->addWidget(magnifierViewer);
    setLayout(vertLayout);
}

void MagnifierWidget::updateView(QRect rect)
{
    QElapsedTimer timer;
    timer.start();

    cv::Mat processed = input.clone();
    cv::Rect cvRect(rect.left(), rect.top(), rect.width(), rect.height());
    cv::Mat roi(processed, cvRect);
    if (equalizeRadio->isChecked())
    {
        std::vector<cv::Mat> chans;
        cv::split(roi, chans);
        for (uint i = 0; i < chans.size(); ++i)
            cv::equalizeHist(chans.at(i), chans[i]);
        cv::merge(chans, roi);
        lastRadio = equalizeRadio;
    }
    else if (contrast1Radio->isChecked())
    {
        cv::Mat gray;
        cv::cvtColor(roi, gray, CV_BGR2GRAY);
        cv::Mat lut = Utility::autoLUT(gray, PERCENTILE);
        cv::LUT(roi, lut, roi);
        lastRadio = contrast1Radio;
    }
    else if (contrast2Radio->isChecked())
    {
        std::vector<cv::Mat> chans;
        cv::split(roi, chans);
        for (uint i = 0; i < chans.size(); ++i)
        {
            cv::Mat lut = Utility::autoLUT(chans.at(i), PERCENTILE);
            cv::LUT(chans.at(i), lut, chans[i]);
        }
        cv::merge(chans, roi);
        lastRadio = contrast2Radio;
    }
    else if (retinexRadio->isChecked())
    {
        /*
        const uchar WHITE_MIN_SAT = 40;
        cv::Ptr<cv::xphoto::GrayworldWB> balancer = cv::xphoto::createGrayworldWB();
        balancer->setSaturationThreshold(WHITE_MIN_SAT / 255.0);
        balancer->balanceWhite(region, balanced);
        return true;
        */
        lastRadio = retinexRadio;
    }
    else
    {
        //if (lastRadio)
            lastRadio->setChecked(true);
        return;
    }
    magnifierViewer->updateProcessed(processed);
    emit messageToShow(tr("[INSPECTION::Magnifier] Processing time = %1 ms").arg(timer.elapsed()));
}

void MagnifierWidget::updateRect()
{
    updateView(magnifierViewer->getSceneRect());
}

