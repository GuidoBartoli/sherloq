#include "originalwidget.h"
#include "singleviewer.h"

OriginalWidget::OriginalWidget(const QString &fileName, const cv::Mat &input, QWidget* parent)
    : ToolWidget(parent)
{
    SingleViewer* imageViewer = new SingleViewer(input, cv::Mat(), QFileInfo(fileName).fileName());
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addWidget(imageViewer);
    setLayout(vertLayout);
}
