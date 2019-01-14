#include "contrastwidget.h"
#include "utility.h"
#include "alglib/interpolation.h"
#include "singleviewer.h"

ContrastWidget::ContrastWidget(const cv::Mat &input, QWidget *parent)
    : ToolWidget(parent)
{
    cv::Mat output;
    contrastEnhancement(input, output, parent);
    SingleViewer* contrastViewer = new SingleViewer(input, output);
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addWidget(contrastViewer);
    setLayout(vertLayout);
}

void ContrastWidget::contrastEnhancement2(const cv::Mat &input, cv::Mat &output, QWidget* parent)
{
    const int    BLOCK_SIZE = std::min(input.rows, input.cols) / 128 > 10 ? 128 : 64;
    const int    CUT_OFF    = 8;
    const double PI         = acos(-1);

    cv::Mat color = Utility::padImage(input, BLOCK_SIZE);
    cv::Mat gray;
    cv::cvtColor(color, gray, CV_BGR2GRAY);
    output = cv::Mat(color.rows / BLOCK_SIZE, color.cols / BLOCK_SIZE, CV_32F);

    /*
    int nimages = 1;
    int chans[] = { 0 };
    cv::Mat mask;
    int dims = 1;
    int histsize[] = { 256 };
    float vranges[] = { 0, 256 };
    const float* ranges[] = { vranges };
    */

    QProgressDialog progress(tr("Detecting Contrast Enhancement..."), QString(), 0, color.rows, parent);
    progress.setWindowModality(Qt::WindowModal);
    //double maxDev = std::numeric_limits<double>::min();
    for (int i = 0; i < color.rows; i += BLOCK_SIZE)
    {
        for (int j = 0; j < color.cols; j += BLOCK_SIZE)
        {
            cv::Mat block(gray, cv::Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            //cv::imwrite((QDir::tempPath() + "/block.png").toStdString(), block);
            //cv::Mat block2 = block.clone();

            cv::Mat hist(1, 256, CV_32F, cv::Scalar(0));
            for (int i = 0; i < block.rows; ++i)
                for (int j = 0; j < block.cols; ++j)
                    ++hist.at<float>(0, block.at<uchar>(i, j));

            //double max = std::numeric_limits<double>::min();
            for (int i = 0; i < hist.cols; ++i)
            {
                if (i <= CUT_OFF)
                    hist.at<float>(0, i) *= (1 - cos((PI*(double)i) / CUT_OFF)) / 2;
                else if (i >= 255 - CUT_OFF)
                    hist.at<float>(0, i) *= (1 + cos((PI*((double)i + CUT_OFF - 255)) / CUT_OFF)) / 2;
                //if (hist.at<float>(0, i) > max)
                //    max = hist.at<float>(0, i);
            }
            //hist /= max;
            cv::blur(hist, hist, cv::Size(3, 1));
            cv::normalize(hist, hist, 0, 1, cv::NORM_MINMAX);
            //cv::imwrite((QDir::tempPath() + "/hist.png").toStdString(), hist);

            /*
            cv::Mat hist;
            cv::calcHist(&block2, nimages, chans, mask, hist, dims, histsize, ranges);
            cv::normalize(hist, hist, 0, 1, cv::NORM_MINMAX);
            */

            double dev = 0;
            for (int k = 2; k <= 253; ++k)
            {
                double h0 = hist.at<float>(0, k - 2);
                double h1 = hist.at<float>(0, k - 1);
                double h2 = hist.at<float>(0, k    );
                double h3 = hist.at<float>(0, k + 1);
                double h4 = hist.at<float>(0, k + 2);
                double h_ = h1 + h3 - (h0 + h4) / 2;
                //double h__ = (h1 + h3) / 2;
                dev += fabs(h2 - h_);
            }
            output.at<float>(i / BLOCK_SIZE, j / BLOCK_SIZE) = dev;
            //if (dev > maxDev)
            //    maxDev = dev;
        }
        progress.setValue(i);
    }
    //output /= maxDev;
    cv::medianBlur(output, output, 3);
    cv::resize(output, output, cv::Size(), BLOCK_SIZE, BLOCK_SIZE, cv::INTER_NEAREST);
    cv::normalize(output, output, 0, 255, cv::NORM_MINMAX);
    output.convertTo(output, CV_8U);
    cv::cvtColor(output, output, CV_GRAY2BGR);
}

void ContrastWidget::contrastEnhancement(const cv::Mat &input, cv::Mat &output, QWidget* parent)
{
    const int BLOCK_SIZE = std::min(input.rows, input.cols) / 128 > 10 ? 128 : 64;

    cv::Mat color = Utility::padImage(input, BLOCK_SIZE);
    cv::Mat gray;
    cv::cvtColor(color, gray, CV_BGR2GRAY);
    output = cv::Mat(color.rows / BLOCK_SIZE, color.cols / BLOCK_SIZE, CV_32F);
    // cv::Mat output1(input0.rows, input0.cols, CV_32F);
    // cv::Mat output2(input0.rows, input0.cols, CV_32F);

    QProgressDialog progress(tr("Detecting Contrast Enhancement..."), QString(), 0, color.rows, parent);
    progress.setWindowModality(Qt::WindowModal);

    for (int i = 0; i < color.rows; i += BLOCK_SIZE)
    {
        for (int j = 0; j < color.cols; j += BLOCK_SIZE)
        {
            cv::Mat color_blk(color, cv::Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            cv::Mat gray_blk(gray, cv::Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            double err = computeHistErr(gray_blk);
            // double spl = computeHistSpl(gray_blk);
            double sim = computeChanSim(color_blk);
            // cv::Mat block1(output1, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            // block1 = err;
            // cv::Mat block2(output2, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            // block2 = sim;
            output.at<float>(i / BLOCK_SIZE, j / BLOCK_SIZE) = err * sim;
        }
        progress.setValue(i);
    }
    cv::medianBlur(output, output, 3);
    cv::resize(output, output, cv::Size(), BLOCK_SIZE, BLOCK_SIZE, cv::INTER_NEAREST);
    output.convertTo(output, CV_8U, 255);
    cv::cvtColor(output, output, CV_GRAY2BGR);

    //TODO: Ritagliare output alle stesse dimensioni dell'originale

    // normalize(output1, output1, 0, 255, NORM_MINMAX);
    // normalize(output2, output2, 0, 255, NORM_MINMAX);
    // output1.convertTo(output1, CV_8U);
    // output2.convertTo(output2, CV_8U);
    // imwrite("err.jpg", output1);
    // imwrite("sim.jpg", output2);
}

double ContrastWidget::computeHistErr(const cv::Mat &input)
{
    const int    CUT_OFF = 8;
    const double MAX_ERR = 0.185;
    const double SPL_POW = 0.5;
    const double PI = acos(-1);

    cv::Mat hist(1, 256, CV_32F, cv::Scalar(0));
    for (int i = 0; i < input.rows; ++i)
        for (int j = 0; j < input.cols; ++j)
            ++hist.at<float>(input.at<uchar>(i, j));

    for (int i = 0; i < hist.cols; ++i)
    {
        if (i <= CUT_OFF)
            hist.at<float>(0, i) *= (1 - cos((PI*(double)i) / CUT_OFF)) / 2;
        else if (i >= 255 - CUT_OFF)
            hist.at<float>(0, i) *= (1 + cos((PI*((double)i+CUT_OFF-255)) / CUT_OFF)) / 2;
    }
    cv::normalize(hist, hist, 0, 1, cv::NORM_MINMAX);

    std::vector<double> x(256), y(256);
    for (int i = 0; i < 256; ++i)
    {
        x[i] = i;
        y[i] = hist.at<float>(i);
    }

    cv::Mat hist_dft = histDFT(hist, true);
    double en = 0, ed = 0, max_diff = std::numeric_limits<double>::min();
    for (int i = 0; i < hist.cols; ++i)
    {
        double w = pow(((double)i - 128) / 128, 2);
        en += fabs(w * hist_dft.at<float>(0, i));
        ed += fabs(hist_dft.at<float>(0, i));

        double v = y.at(i);
        x.erase(x.begin() + i);
        y.erase(y.begin() + i);
        alglib::spline1dinterpolant spline;
        alglib::real_1d_array ax, ay;
        ax.setcontent(x.size(), &(x[0]));
        ay.setcontent(y.size(), &(y[0]));
        alglib::spline1dbuildcubic(ax, ay, spline);
        double s = alglib::spline1dcalc(spline, i);
        double diff = fabs(v - s);
        if (diff > max_diff)
            max_diff = diff;
        x.insert(x.begin() + i, i);
        y.insert(y.begin() + i, v);
    }

    double err = en / ed;
    if (err > MAX_ERR)
        err = 1;
    else
        err /= MAX_ERR;

    return (err * pow(max_diff, SPL_POW));
}

cv::Mat ContrastWidget::histDFT(const cv::Mat &hist, bool full)
{
    int len = hist.cols;
    cv::Mat complex, planes[] = {cv::Mat_<float>(hist), cv::Mat::zeros(1, len, CV_32F)};
    cv::merge(planes, 2, complex);
    cv::dft(complex, complex);
    cv::split(complex, planes);
    cv::magnitude(planes[0], planes[1], planes[0]);
    cv::Mat dft_mag;
    if (full)
        cv::hconcat(planes[0].colRange(len/2, len), planes[0].colRange(0, len/2), dft_mag);
    else
        dft_mag = planes[0].colRange(0, len/2).clone();
    cv::normalize(dft_mag, dft_mag, 0, 1, cv::NORM_MINMAX);
    return dft_mag;
}

double ContrastWidget::computeChanSim(const cv::Mat &input)
{
    const double MAX_SIM = 0.75;

    std::vector<cv::Mat> bgr, deriv;
    cv::split(input, bgr);
    cv::Mat kx, ky;
    getDerivKernels(kx, ky, 1, 1, 1);
    for (uint c = 0; c < bgr.size(); ++c)
    {
        cv::Mat flt;
        sepFilter2D(bgr.at(c), flt, CV_32F, kx, ky);
        deriv.push_back(flt);
    }

    double s = 0, m = 0;
    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            float b = deriv.at(0).at<float>(i, j);
            float g = deriv.at(1).at<float>(i, j);
            float r = deriv.at(2).at<float>(i, j);
            s += (fabs(g - r) + fabs(g - b) + fabs(r - b));
            m += (fabs(r) + fabs(g) + fabs(b));
        }
    }
    s /= (3 * input.rows * input.cols);
    m /= (3 * input.rows * input.cols);
    double sim = s / m;

    return (sim > MAX_SIM ? 1 : sim / MAX_SIM);
}

