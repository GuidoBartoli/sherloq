#include "ghostwidget.h"
#include "utility.h"

GhostWidget::GhostWidget(const cv::Mat &input, QWidget* parent)
    : ToolWidget(parent)
{
    cv::Mat qualities, losses;
    int n, l1, l2, q1, q2;
    float p1, p2;
    computeGhosts(input, qualities, losses, n, l1, p1, l2, p2, q1, q2, parent);
    computeHist(input, histogram, spectrum);

    QLineSeries* lossSeries = new QLineSeries();
    for (int i = 0; i < qualities.cols; ++i)
        lossSeries->append(qualities.at<uchar>(0, i), losses.at<float>(0, i) * 100);
    QChart *lossChart = new QChart();
    lossChart->legend()->hide();
    lossChart->setTitle(tr("Residual Error -vs- Compression Level"));
    lossChart->addSeries(lossSeries);
    lossChart->createDefaultAxes();
    lossChart->axisX()->setRange(0, 100);
    lossChart->axisX()->setTitleText(tr("quality (%)"));
    static_cast<QValueAxis*>(lossChart->axisX())->setTickCount(11);
    static_cast<QValueAxis*>(lossChart->axisX())->setLabelFormat("%d");
    lossChart->axisY()->setTitleText(tr("loss (%)"));
    lossChart->setMinimumSize(600, 400);
    QFont font = lossChart->titleFont();
    font.setBold(true);
    lossChart->setTitleFont(font);
    QChartView* lossView = new QChartView(lossChart);
    lossView->setRenderHint(QPainter::Antialiasing);

    QLabel* resultLabel = new QLabel();
    if (l1 < 0)
        resultLabel->setText(tr("Uncompressed (P = %1%, QF = [%2%;%3%])").arg(p1).arg(q1).arg(q2));
    else
    {
        if (l2 < 0)
            resultLabel->setText(tr("Single compression (QF = %1%, P = %2%)").arg(l1).arg(p1));
        else
            resultLabel->setText(tr("Multiple compression [N%1] (QF1 = %2%, P1 = %3%) ==> (QF2 = %4%, P2 = %5%)").arg(n).arg(l1).arg(p1).arg(l2).arg(p2));
    }
    font = resultLabel->font();
    font.setBold(true);
    resultLabel->setFont(font);

    QTableWidget* dctTable = new QTableWidget(8, 8);
    cv::Mat color(1, 1, CV_8UC3);
    for (int k = 0; k < TABLE_SIZE; ++k)
    {
        int i = Utility::ZIG_ZAG.at<int>(k, 0);
        int j = Utility::ZIG_ZAG.at<int>(k, 1);
        QTableWidgetItem* item = new QTableWidgetItem();
        item->setText(QString::number(k + 1));
        item->setTextAlignment(Qt::AlignRight);
        color.at<cv::Vec3b>(0, 0) = cv::Vec3b((uchar)((k / 64.0) * 160.0) + 10, 255, 255);
        cv::cvtColor(color, color, CV_HSV2BGR);
        cv::Vec3b vec = color.at<cv::Vec3b>(0, 0);
        item->setBackgroundColor(QColor(vec[0], vec[1], vec[2]));
        dctTable->setItem(i, j, item);
    }
    //dctTable->resizeColumnsToContents();
    dctTable->resizeRowsToContents();
    dctTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    dctTable->setSelectionMode(QAbstractItemView::SingleSelection);
    for (int i = 0; i < dctTable->columnCount(); ++i)
        dctTable->setColumnWidth(i, 30);
    //dctTable->horizontalHeader()->hide();
    //dctTable->verticalHeader()->hide();
    connect(dctTable, &QTableWidget::itemClicked,
            this,     &GhostWidget::changeMode);

    dctChart = new QChart();
    QChartView* dctView = new QChartView(dctChart);
    dctView->setRenderHint(QPainter::Antialiasing);
    dftChart = new QChart();
    QChartView* dftView = new QChartView(dftChart);
    dftView->setRenderHint(QPainter::Antialiasing);
    changeMode(new QTableWidgetItem("1"));

    QGridLayout* gridLayout = new QGridLayout();
    gridLayout->addWidget(lossView, 0, 0, 1, 2);
    gridLayout->addWidget(resultLabel, 1, 0, 1, 2);
    gridLayout->addWidget(dctTable, 2, 0);
    gridLayout->addWidget(dctView, 2, 1);
    //gridLayout->addWidget(dftView, 2, 2);
    setLayout(gridLayout);

}

void GhostWidget::computeGhosts(const cv::Mat &image, cv::Mat &qualities,
                                cv::Mat &losses, int &n, int &l1, float &p1,
                                int &l2, float &p2, int &q1, int &q2, QWidget* parent)
{
    cv::Mat original;
    cv::cvtColor(image, original, CV_BGR2GRAY);
    original.convertTo(original, CV_32F);

    qualities.create(1, (MAX_QUALITY - MIN_QUALITY) / QUALITY_STEP + 1, CV_8U);
    int i = 0;
    for (int q = MAX_QUALITY; q >= MIN_QUALITY; q -= QUALITY_STEP)
        qualities.at<uchar>(0, i++) = q;
    losses.create(1, qualities.cols, CV_32F);

    QProgressDialog progress(tr("Computing residuals..."), QString(), 0, qualities.cols, parent);
    progress.setWindowModality(Qt::WindowModal);
    for (int i = 0; i < qualities.cols; ++i)
    {
        progress.setValue(i);
        cv::Mat compressed = Utility::jpegCompress(original, qualities.at<uchar>(0, i), false);
        cv::imwrite("/tmp/comp.jpg", compressed);
        compressed.convertTo(compressed, CV_32F);
        cv::Scalar loss = cv::mean(cv::abs(original - compressed));
        losses.at<float>(0, i) = loss[0];
        //qDebug() << qualities.at<uchar>(0, i) << loss[0];
    }
    cv::normalize(losses, losses, 0, 1, cv::NORM_MINMAX);

    cv::Mat response;
    cv::Mat kernel = (cv::Mat_<float>(1, KERNEL_SIZE) << 0, 0, -1, 0, 0);
    cv::matchTemplate(losses, kernel, response, CV_TM_CCOEFF_NORMED);
    cv::Mat padding(1, KERNEL_SIZE / 2, CV_32F, cv::Scalar(0));
    cv::hconcat(padding, response, response);
    cv::hconcat(response, padding, response);
    double peak_val;
    cv::Point peak_loc;
    cv::minMaxLoc(response, 0, &peak_val, 0, &peak_loc);

    if (peak_val < THR_PEAK)
    {
        l1 = l2 = p2 = -1;
        p1 = (peak_val - THR_PEAK) / (MIN_PEAK - THR_PEAK) * 100;
        n = 0;
    }
    else
    {
        l1 = qualities.at<uchar>(peak_loc.x);
        p1 = (peak_val - THR_PEAK) / (MAX_PEAK - THR_PEAK) * 100;
        response.at<float>(0, peak_loc.x) = response.at<float>(0, peak_loc.x - 1) =
                response.at<float>(0, peak_loc.x + 1) = -1;
        cv::minMaxLoc(response, 0, &peak_val, 0, &peak_loc);
        if (peak_val > THR_PEAK)
        {
            l2 = qualities.at<uchar>(peak_loc.x);
            p2 = (peak_val - THR_PEAK) / (MAX_PEAK - THR_PEAK) * 100;
            if (p1 > p2)
            {
                std::swap(l1, l2);
                std::swap(p1, p2);
            }
            n = 2;
            for (int i = 0; i < response.cols; ++i)
                if (response.at<float>(0, i) > THR_PEAK)
                    ++n;
        }
        else
        {
            l2 = p2 = -1;
            n = 1;
        }
    }
    if (p1 > 100)
        p1 = 100;
    if (p2 > 100)
        p2 = 100;
    q1 = MIN_QUALITY + KERNEL_SIZE / 2;
    q2 = MAX_QUALITY - KERNEL_SIZE / 2;
}

void GhostWidget::computeHist(const cv::Mat &image, cv::Mat &histogram, cv::Mat &spectrum)
{
    cv::Mat gray;
    cv::cvtColor(image, gray, CV_BGR2GRAY);
    cv::Mat padded = Utility::padImage(gray, DCT_SIZE);
    padded.convertTo(padded, CV_32F, 1, -CHAR_SHIFT);

    histogram.create(TABLE_SIZE, 2*MAX_COEFF, CV_32F);
    histogram.setTo(0);
    cv::Mat coeffs(DCT_SIZE, DCT_SIZE, CV_32F);
    for (int i = 0; i < padded.rows - DCT_SIZE; i += DCT_SIZE)
    {
        for (int j = 0; j < padded.cols - DCT_SIZE; j += DCT_SIZE)
        {
            cv::Mat block(padded, cv::Rect(j, i, DCT_SIZE, DCT_SIZE));
            cv::dct(block, coeffs);
            for (int k = 0; k < TABLE_SIZE; ++k)
            {
                int i0 = Utility::ZIG_ZAG.at<int>(k, 0);
                int j0 = Utility::ZIG_ZAG.at<int>(k, 1);
                float coeff = coeffs.at<float>(i0, j0);
                if (coeff < -MAX_COEFF || coeff >= +MAX_COEFF)
                    continue;
                ++histogram.at<float>(k, coeff + MAX_COEFF);
            }
        }
    }
    cv::normalize(histogram, histogram, 0, 1, cv::NORM_MINMAX);

    spectrum.create(histogram.rows, histogram.cols, CV_32F);
    for (int i = 0; i < histogram.rows; ++i)
    {
        cv::Mat row = histogram.row(i).clone();

        /*
        cv::Mat row2;                                   //expand input image to optimal size
        int m = cv::getOptimalDFTSize( row.rows );
        int n = cv::getOptimalDFTSize( row.cols );      // on the border add zero pixels
        copyMakeBorder(row, row2, 0, m - row.rows, 0, n - row.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));
        */

        cv::Mat planes[] = {cv::Mat_<float>(row), cv::Mat::zeros(row.size(), CV_32F)};
        cv::Mat complexI;
        cv::merge(planes, 2, complexI);                 // Add to the expanded another plane with zeros
        cv::dft(complexI, complexI);                    // this way the result may fit in the source matrix
        cv::split(complexI, planes);                    // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
        cv::magnitude(planes[0], planes[1], planes[0]); // planes[0] = magnitude
        cv::Mat magI = planes[0];
        magI += cv::Scalar::all(1);                     // switch to logarithmic scale
        cv::log(magI, magI);
        //magI = magI(cv::Rect(0, 0, magI.cols & -2, magI.rows & -2));
        cv::normalize(magI, magI, 0, 1, cv::NORM_MINMAX);
        spectrum.row(i) = magI;
    }
    //Utility::saveImage("spectrum", spectrum);
}

void GhostWidget::changeMode(QTableWidgetItem* item)
{
    dctChart->removeAllSeries();
    QLineSeries* dctSeries = new QLineSeries();
    int k = item->text().toInt() - 1;
    int i0 = Utility::ZIG_ZAG.at<int>(k, 0) + 1;
    int j0 = Utility::ZIG_ZAG.at<int>(k, 1) + 1;
    for (int i = 0; i < 2*MAX_COEFF; ++i)
        dctSeries->append(i - MAX_COEFF, histogram.at<float>(k, i) * 100);

    dctChart->addSeries(dctSeries);
    dctChart->createDefaultAxes();
    dctChart->legend()->hide();
    dctChart->axisX()->setRange(-MAX_COEFF, +MAX_COEFF);
    dctChart->axisX()->setTitleText(tr("value"));
    dctChart->axisY()->setTitleText(tr("amount"));
    dctChart->setTitle(tr("DCT Histogram (%1,%2)").arg(i0).arg(j0));

    dftChart->removeAllSeries();
    QLineSeries* dftSeries = new QLineSeries();
    for (int i = 0; i < 2*MAX_COEFF; ++i)
        dftSeries->append(i - MAX_COEFF, spectrum.at<float>(k, i));
    dftChart->addSeries(dftSeries);
    dftChart->createDefaultAxes();
    dftChart->legend()->hide();
    dftChart->axisX()->setRange(-MAX_COEFF, +MAX_COEFF);
    dftChart->axisX()->setTitleText(tr("value"));
    dftChart->axisY()->setTitleText(tr("amount"));
    dftChart->setTitle(tr("DFT Histogram (%1,%2)").arg(i0).arg(j0));
}


