#include "minmaxwidget.h"
#include "singleviewer.h"

MinMaxWidget::MinMaxWidget(const cv::Mat &input, QWidget *parent)
    : ToolWidget(parent)
{
    cv::Mat output = cv::Mat::zeros(input.rows, input.cols, CV_8UC3);
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
                output.at<cv::Vec3b>(i, j) = cv::Vec3b(0, 0, 255);
            else if (norm0 < min_norm)
                output.at<cv::Vec3b>(i, j) = cv::Vec3b(0, 255, 0);
        }
    }
    SingleViewer* minMaxViewer = new SingleViewer(input, output);
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addWidget(minMaxViewer);
    setLayout(vertLayout);
}
