#include "utility.h"

/*
QImage Utility::cv2img(const cv::Mat &input)
{
    return QImage(input.data, input.cols, input.rows, input.step, QImage::Format_RGB888).rgbSwapped();
}

cv::Mat Utility::img2cv(const QImage &input)
{
    return cv::Mat(input.height(), input.width(), CV_8UC3, input.rgbSwapped().scanLine(0));
}
*/

cv::Mat Utility::padImage(const cv::Mat &image, int block_size)
{
    int right = 0, bottom = 0;
    if (image.cols % block_size != 0)
        right = ((image.cols / block_size) * block_size + block_size) - image.cols;
    if (image.rows % block_size != 0)
        bottom = ((image.rows / block_size) * block_size + block_size) - image.rows;
    if (right == 0 && bottom == 0)
        return image.clone();
    cv::Mat padded;
    if (image.channels() == 1)
        cv::copyMakeBorder(image, padded, 0, bottom, 0, right, cv::BORDER_CONSTANT, cv::Scalar(0));
    else
        cv::copyMakeBorder(image, padded, 0, bottom, 0, right, cv::BORDER_CONSTANT, cv::Scalar(0,0,0));
    return padded;
}

/*
QString Utility::formatRange(int val)
{
    QString string;
    if (val > 0)
        string.append("+");
    else if (val < 0)
        string.append("-");
    else
        string.append(" ");
    if (val < 10)
        string.append("  ");
    else if (val < 100)
        string.append(" ");
    string.append(QString::number(abs(val)));
    return string;
}
*/

cv::Mat Utility::jpegCompress(const cv::Mat &input, int quality, bool color)
{
    std::vector<int> params;
    params.push_back(CV_IMWRITE_JPEG_QUALITY);
    params.push_back(quality);
    std::vector<uchar> buffer;
    cv::imencode(".jpg", input, buffer, params);
    return cv::imdecode(buffer, color ? CV_LOAD_IMAGE_COLOR : CV_LOAD_IMAGE_GRAYSCALE);
}

/*
cv::Mat Utility::buildContrastLUT(uchar low, uchar high)
{
    cv::Mat lut(256, 1, CV_8U);
    if (low > high)
        std::swap(low, high);
    for (int i = 0; i < 256; ++i)
    {
        if (i <= low)
            lut.at<uchar>(i, 0) = 0;
        else if (i >= high)
            lut.at<uchar>(i, 0) = 255;
        else
            lut.at<uchar>(i, 0) = (uchar)((i - low) * (255.0 / (high - low)));
    }
    return lut;
}
*/

cv::Mat Utility::createLUT(int low, int high)
{
    if (abs(low) + abs(high) >= 255)
        low = high = 127;

    double x1, y1, x2, y2;
    if (low >= 0)
    {
        x1 = low;
        y1 = 0;
    }
    else
    {
        x1 = 0;
        y1 = -low;
    }
    if (high >= 0)
    {
        x2 = 255 - high;
        y2 = 255;
    }
    else
    {
        x2 = 255;
        y2 = 255 + high;
    }

    cv::Mat lut(256, 1, CV_8U);
    for (int i = 0; i <= 255; ++i)
    {
        if (low >= 0 && i <= x1)
            lut.at<uchar>(i, 0) = 0;
        else if (high >= 0 && i >= x2)
            lut.at<uchar>(i, 0) = 255;
        else
            lut.at<uchar>(i, 0) = (uchar)((i*y1-x2*y1-i*y2+x1*y2)/(x1-x2));
        //qDebug() << i << lut.at<uchar>(i, 0);
    }
    return lut;
}

cv::Mat Utility::autoLUT(const cv::Mat &channel, float percentile)
{
    cv::Mat hist = getHist(channel);

    float countThr = channel.total() * percentile;
    int lowCount = 0;
    int highCount = 0;
    uchar low = 0;
    uchar high = 255;
    for (int i = 0; i < 256; ++i)
    {
        lowCount += hist.at<float>(i, 0);
        if (lowCount >= countThr)
        {
            low = i;
            break;
        }
    }
    for (int i = 255; i >= 0; --i)
    {
        highCount += hist.at<float>(i, 0);
        if (highCount >= countThr)
        {
            high = i;
            break;
        }
    }

    return Utility::createLUT(low, 255 - high);
}

void Utility::changeFont(QWidget* widget, bool bold, bool italic)
{
    QFont font = widget->font();
    font.setBold(bold);
    font.setItalic(italic);
    widget->setFont(font);
}

void Utility::saveImage(QString name, const cv::Mat &image)
{
    cv::Mat temp;
    cv::normalize(image, temp, 0, 255, cv::NORM_MINMAX);
    temp.convertTo(temp, CV_8U);
    cv::imwrite((QDir::tempPath() + QString("/") + name +
                 QString(".jpg")).toStdString(), temp);
}

cv::Mat Utility::getHist(const cv::Mat &channel, bool normalize)
{
    /*
    int nimages = 1;
    int chans[] = { 0 };
    cv::Mat mask;
    cv::Mat hist;
    int dims = 1;
    int histsize[] = { 256 };
    float vranges[] = { 0, 256 };
    const float* ranges[] = { vranges };
    cv::calcHist(&channel, nimages, chans, mask, hist, dims, histsize, ranges);
    if (normalize)
        cv::normalize(hist, hist, 0, 1, cv::NORM_MINMAX);
        //hist /= channel.total();
    */
    cv::Mat hist(256, 1, CV_32F, cv::Scalar(0));
    for (int i = 0; i < channel.rows; ++i)
        for (int j = 0; j < channel.cols; ++j)
            ++hist.at<float>(channel.at<uchar>(i, j), 0);
    if (normalize)
        cv::normalize(hist, hist, 0, 1, cv::NORM_MINMAX);

    return hist;
}
