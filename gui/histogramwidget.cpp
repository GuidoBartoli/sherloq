#include "histogramwidget.h"
#include "utility.h"

HistogramWidget::HistogramWidget(const cv::Mat &image, QWidget *parent)
    : ToolWidget(parent)
{
    cv::Mat gray;
    cv::cvtColor(image, gray, CV_BGR2GRAY);
    std::vector<cv::Mat> bgr;
    cv::split(image, bgr);
    cv::Mat red, green, blue;
    red = bgr.at(2);
    green = bgr.at(1);
    blue = bgr.at(0);

    cv::Mat lumHist = Utility::getHist(gray, true);
    cv::Mat redHist = Utility::getHist(red, true);
    cv::Mat greenHist = Utility::getHist(green, true);
    cv::Mat blueHist = Utility::getHist(blue, true);

    QLineSeries* lumLine   = new QLineSeries();
    QLineSeries* redLine   = new QLineSeries();
    QLineSeries* greenLine = new QLineSeries();
    QLineSeries* blueLine  = new QLineSeries();
    QLineSeries* baseLine  = new QLineSeries();
    for (int i = 0; i < 256; ++i)
    {
        lumLine->append(i, lumHist.at<float>(i, 0) * 100);
        redLine->append(i, redHist.at<float>(i, 0) * 100);
        qDebug() << (int)(redHist.at<float>(i, 0) * 100);
        greenLine->append(i, greenHist.at<float>(i, 0) * 100);
        blueLine->append(i, blueHist.at<float>(i, 0) * 100);
        baseLine->append(i, 0);
    }

    uchar alpha = 128;
    lumArea = new QAreaSeries(lumLine, baseLine);
    QColor color = Qt::lightGray;
    color.setAlpha(alpha);
    lumArea->setBrush(color);
    lumArea->setPen(QPen(Qt::white));

    redArea = new QAreaSeries(redLine, baseLine);
    color = Qt::darkRed;
    color.setAlpha(alpha);
    redArea->setBrush(color);
    redArea->setPen(QPen(Qt::red));

    greenArea = new QAreaSeries(greenLine, baseLine);
    color = Qt::darkGreen;
    color.setAlpha(alpha);
    greenArea->setBrush(color);
    greenArea->setPen(QPen(Qt::green));

    blueArea = new QAreaSeries(blueLine, baseLine);
    color = Qt::darkBlue;
    color.setAlpha(alpha);
    blueArea->setBrush(color);
    blueArea->setPen(QPen(Qt::blue));

    histChart = new QChart();
    histChart->addSeries(lumArea);
    histChart->createDefaultAxes();
    histChart->setPlotAreaBackgroundBrush(Qt::darkGray);
    histChart->setPlotAreaBackgroundVisible(true);
    histChart->setBackgroundVisible(false);
    histChart->legend()->hide();
    histChart->axisX()->setRange(0, 255);
    histChart->axisX()->setTitleText(tr("level (0-255)"));
    static_cast<QValueAxis*>(histChart->axisX())->setTickCount(5);
    static_cast<QValueAxis*>(histChart->axisX())->setMinorTickCount(1);
    static_cast<QValueAxis*>(histChart->axisX())->setLabelFormat("%d");
    histChart->axisY()->setRange(0, 100);
    histChart->axisY()->setTitleText(tr("count (%)"));

    QChartView* histView = new QChartView(histChart);
    histView->setRenderHint(QPainter::Antialiasing);

    QRadioButton* lumRadio   = new QRadioButton(tr("Luminance"));
    QRadioButton* redRadio   = new QRadioButton(tr("Red"));
    QRadioButton* greenRadio = new QRadioButton(tr("Green"));
    QRadioButton* blueRadio  = new QRadioButton(tr("Blue"));
    QRadioButton* compRadio  = new QRadioButton(tr("Composite"));
    QCheckBox*    logCheck   = new QCheckBox(tr("Logarithmic"));

    QHBoxLayout* horizLayout = new QHBoxLayout();
    horizLayout->addWidget(new QLabel(tr("Channel:")));
    horizLayout->addWidget(lumRadio);
    horizLayout->addWidget(redRadio);
    horizLayout->addWidget(greenRadio);
    horizLayout->addWidget(blueRadio);
    horizLayout->addWidget(compRadio);
    horizLayout->addWidget(logCheck);
    horizLayout->addStretch();

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(horizLayout);
    vertLayout->addWidget(histView);
    setLayout(vertLayout);
    setMinimumHeight(500);

    connect(lumRadio,   &QRadioButton::toggled,
            this,       &HistogramWidget::setLumHist);
    connect(redRadio,   &QRadioButton::toggled,
            this,       &HistogramWidget::setRedHist);
    connect(greenRadio, &QRadioButton::toggled,
            this,       &HistogramWidget::setGreenHist);
    connect(blueRadio,  &QRadioButton::toggled,
            this,       &HistogramWidget::setBlueHist);
    connect(compRadio,  &QRadioButton::toggled,
            this,       &HistogramWidget::setCompHist);
    connect(logCheck,   &QCheckBox::toggled,
            this,       &HistogramWidget::setLogHist);

    logCheck->setEnabled(false);
    lumRadio->toggle();
}

void HistogramWidget::setLumHist(bool toggled)
{
    if (toggled)
    {
        QList<QAbstractSeries*> series = histChart->series();
        for (int i = 0; i < series.length(); ++i)
            histChart->removeSeries(series.at(i));
        histChart->addSeries(lumArea);
    }
}

void HistogramWidget::setRedHist(bool toggled)
{
    if (toggled)
    {
        QList<QAbstractSeries*> series = histChart->series();
        for (int i = 0; i < series.length(); ++i)
            histChart->removeSeries(series.at(i));
        histChart->addSeries(redArea);
    }
}

void HistogramWidget::setGreenHist(bool toggled)
{
    if (toggled)
    {
        QList<QAbstractSeries*> series = histChart->series();
        for (int i = 0; i < series.length(); ++i)
            histChart->removeSeries(series.at(i));
        histChart->addSeries(greenArea);
    }
}

void HistogramWidget::setBlueHist(bool toggled)
{
    if (toggled)
    {
        QList<QAbstractSeries*> series = histChart->series();
        for (int i = 0; i < series.length(); ++i)
            histChart->removeSeries(series.at(i));
        histChart->addSeries(blueArea);
    }
}

void HistogramWidget::setCompHist(bool toggled)
{
    if (toggled)
    {
        QList<QAbstractSeries*> series = histChart->series();
        for (int i = 0; i < series.length(); ++i)
            histChart->removeSeries(series.at(i));
        histChart->addSeries(redArea);
        histChart->addSeries(greenArea);
        histChart->addSeries(blueArea);
    }
}

void HistogramWidget::setLogHist(bool checked)
{
    /*
    if (checked)
    {
        histChart->removeAxis(linearAxis);
        histChart->addAxis(logAxis, Qt::AlignLeft);
        histChart->series().at(0)->attachAxis(logAxis);
    }
    else
    {
        histChart->removeAxis(logAxis);
        histChart->addAxis(linearAxis, Qt::AlignLeft);
        histChart->series().at(0)->attachAxis(linearAxis);
    }
    */
}
