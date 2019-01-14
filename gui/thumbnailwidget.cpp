#include "thumbnailwidget.h"
#include <opencv2/opencv.hpp>
#include "singleviewer.h"

ThumbnailWidget::ThumbnailWidget(const QString &fileName, QWidget* parent)
    : ToolWidget(parent)
{
    ExifTool* exif = new ExifTool();
    //TagInfo* info = exif->ImageInfo(fileName.toStdString().c_str(), "-b\n-PreviewImage");
    TagInfo* info = exif->ImageInfo(fileName.toStdString().c_str(), "-b\n-ThumbnailImage");
    if (!info)
    {
        showError(0, exif, info);
        return;
    }
    if (!info->next)
    {
        //delete info;
        //info = exif->ImageInfo(fileName.toStdString().c_str(), "-b\n-ThumbnailImage");
        showError(1, exif, info);
        return;
    }
    char* data = info->next->value;
    int len = info->next->valueLen;
    std::vector<char> buffer(data, data + len);
    cv::Mat thumb = cv::imdecode(buffer, CV_LOAD_IMAGE_COLOR);
    if (thumb.empty())
    {
        showError(2, exif, info);
        return;
    }

    cv::Mat image = cv::imread(fileName.toStdString(), CV_LOAD_IMAGE_COLOR);
    //cv::cvtColor(thumb, thumb, CV_BGR2GRAY);
    double imageRatio = (double)image.cols / image.rows;
    double thumbRatio = (double)thumb.cols / thumb.rows;
    cv::Size size;
    if (thumbRatio < imageRatio)
    {
        size.width = thumb.cols;
        size.height = size.width / imageRatio;
    }
    else if (thumbRatio > imageRatio)
    {
        size.height = thumb.rows;
        size.width = size.height * imageRatio;
    }
    else
    {
        size.width = thumb.cols;
        size.height = thumb.rows;
    }

    cv::Mat filtered, resized, padded;
    cv::GaussianBlur(image, filtered, cv::Size(7, 7), 0);
    cv::resize(filtered, resized, size, 0, 0, cv::INTER_AREA);
    int top = (thumb.rows - resized.rows) / 2;
    int bottom = top;
    int left = (thumb.cols - resized.cols) / 2;
    int right = left;
    cv::copyMakeBorder(resized, padded, top, bottom, left, right, cv::BORDER_CONSTANT);
    if (padded.size != thumb.size)
        cv::resize(padded, padded, cv::Size(thumb.cols, thumb.rows));
//    cv::imwrite(QDir::tempPath().toStdString() + "/resized.jpg", resized);
//    cv::imwrite(QDir::tempPath().toStdString() + "/thumb.jpg", thumb);
//    cv::imwrite(QDir::tempPath().toStdString() + "/padded.jpg", padded);
    cv::Mat diff = cv::abs(thumb - padded);
    cv::cvtColor(diff, diff, CV_BGR2GRAY);

    SingleViewer* thumbView = new SingleViewer(thumb, diff);
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addWidget(thumbView);
    setLayout(vertLayout);
    delete info;
    delete exif;
}

void ThumbnailWidget::showError(int type, ExifTool* exif, TagInfo* info)
{
    QLabel* errorLabel = new QLabel();
    if (type == 0)
        errorLabel->setText(tr("Error while extracting information from file!"));
    else if (type == 1)
        errorLabel->setText(tr("File does not contain any thumbnail image"));
    else if (type == 2)
        errorLabel->setText(tr("Error while reading thumbnail image!"));
    else
        errorLabel->setText(tr("Unknown error!"));
    errorLabel->setAlignment(Qt::AlignCenter);
    QFont font = errorLabel->font();
    font.setItalic(true);
    errorLabel->setFont(font);
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addWidget(errorLabel);
    setLayout(vertLayout);
    delete info;
    delete exif;
}
