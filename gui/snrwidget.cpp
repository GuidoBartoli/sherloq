#include "snrwidget.h"
#include "singleviewer.h"

SnrWidget::SnrWidget(const cv::Mat &input, QWidget* parent )
    : ToolWidget(parent)
{
    const int BLOCK_SIZE = 8;

    cv::Mat output;
    cv::Mat min_max;
    minMaxDeviation(input, min_max);
    output.create(input.rows, input.cols, CV_32F);
    output.setTo(255);
    cv::Mat mask(input.rows, input.cols, CV_8U, cv::Scalar(0));
    for (int i = 0; i < min_max.rows - BLOCK_SIZE; i += BLOCK_SIZE)
    {
        for (int j = 0; j < min_max.cols - BLOCK_SIZE; j += BLOCK_SIZE)
        {
            cv::Mat block1(min_max, cv::Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            cv::Scalar avg, stddev;
            cv::meanStdDev(block1, avg, stddev);
            cv::Mat block2(output, cv::Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            block2 = stddev[1] + stddev[2];
            cv::Mat block3(mask, cv::Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            block3 = 1;
        }
    }
    cv::normalize(output, output, 0, 255, cv::NORM_MINMAX, -1, mask);
    output = 255 - output;
    output.convertTo(output, CV_8U);
    cv::cvtColor(output, output, CV_GRAY2BGR);

    SingleViewer* snrViewer = new SingleViewer(input, output);
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addWidget(snrViewer);
    setLayout(vertLayout);
}

void SnrWidget::minMaxDeviation(const cv::Mat &input, cv::Mat &output)
{
    const cv::Vec3b MAX_COLOR(0, 0, 255);
    const cv::Vec3b MIN_COLOR(0, 255, 0);

    output = cv::Mat::zeros(input.rows, input.cols, CV_8UC3);
    for (int i = 1; i < input.rows - 1; ++i)
    {
        for (int j = 1; j < input.cols - 1; ++j)
        {
            float norm0 = cv::norm(input.at<cv::Vec3b>(i, j));
            float norm1 = cv::norm(input.at<cv::Vec3b>(i, j - 1));
            float norm2 = cv::norm(input.at<cv::Vec3b>(i - 1, j));
            float norm3 = cv::norm(input.at<cv::Vec3b>(i, j + 1));
            float norm4 = cv::norm(input.at<cv::Vec3b>(i + 1, j));
            float max_norm = std::max(std::max(std::max(norm1, norm2), norm3), norm4);
            float min_norm = std::min(std::min(std::min(norm1, norm2), norm3), norm4);
            if (norm0 > max_norm)
                output.at<cv::Vec3b>(i, j) = MAX_COLOR;
            else if (norm0 < min_norm)
                output.at<cv::Vec3b>(i, j) = MIN_COLOR;
        }
    }
}
