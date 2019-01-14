#ifndef UTILITY_H
#define UTILITY_H

#include <opencv2/opencv.hpp>
#include <QtWidgets>

namespace Utility
{
    const int DCT_SIZE    = 8;
    const int TABLE_SIZE  = DCT_SIZE * DCT_SIZE;
    const cv::Mat ZIG_ZAG = (cv::Mat_<int>(TABLE_SIZE, 2) <<
        0, 0,   0, 1,   1, 0,   2, 0,   1, 1,   0, 2,   0, 3,   1, 2,
        2, 1,   3, 0,   4, 0,   3, 1,   2, 2,   1, 3,   0, 4,   0, 5,
        1, 4,   2, 3,   3, 2,   4, 1,   5, 0,   6, 0,   5, 1,   4, 2,
        3, 3,   2, 4,   1, 5,   0, 6,   0, 7,   1, 6,   2, 5,   3, 4,
        4, 3,   5, 2,   6, 1,   7, 0,   7, 1,   6, 2,   5, 3,   4, 4,
        3, 5,   2, 6,   1, 7,   2, 7,   3, 6,   4, 5,   5, 4,   6, 3,
        7, 2,   7, 3,   6, 4,   5, 5,   4, 6,   3, 7,   4, 7,   5, 6,
        6, 5,   7, 4,   7, 5,   6, 6,   5, 7,   6, 7,   7, 6,   7, 7);

    //QImage cv2img(const cv::Mat &input);
    //cv::Mat img2cv(const QImage &input);
    cv::Mat jpegCompress(const cv::Mat &input, int quality, bool color = true);
    //cv::Mat buildContrastLUT(uchar low = 0, uchar high = 90);
    cv::Mat createLUT(int low = 0, int high = 165);
    cv::Mat autoLUT(const cv::Mat &channel, float percentile);
    //QString formatRange(int val);
    cv::Mat padImage(const cv::Mat &image, int block_size);
    void changeFont(QWidget* widget, bool bold = false, bool italic = false);
    void saveImage(QString name, const cv::Mat &image);
    cv::Mat getHist(const cv::Mat &channel, bool normalize = false);
}

#endif // UTILITY_H
