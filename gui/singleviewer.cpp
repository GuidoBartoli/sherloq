#include "singleviewer.h"
#include "cvmatandqimage.h"
#include "utility.h"

SingleViewer::SingleViewer(const cv::Mat &original, const cv::Mat &processed, const QString &title, QWidget* parent)
    : QWidget(parent),
      original(original),
      processed(processed)
{
    QGraphicsScene* imageScene = new QGraphicsScene(this);
    imageItem = new QGraphicsPixmapItem(QPixmap::fromImage(QtOcv::mat2Image(original)));
    imageScene->setBackgroundBrush(Qt::darkGray);
    imageScene->addItem(imageItem);
    imageView = new DynamicView(imageScene);

    QLabel* viewLabel = new QLabel(tr("View:"));
    originalRadio = new QRadioButton(tr("Original"));
    originalRadio->setToolTip(tr("Show the reference image for comparison"));
    processRadio = new QRadioButton(tr("Processed"));
    processRadio->setToolTip(tr("Show the image processed with the current tool"));
    zoomLabel = new QLabel();
    QToolButton* fullButton = new QToolButton();
    fullButton->setText(tr("Full"));
    QToolButton* fitButton = new QToolButton();
    fitButton->setText(tr("Fit"));
    //double mpxs = original.rows * original.cols / 1.0e6;
    int width = original.cols;
    int height = original.rows;
    //QLabel* sizeLabel = new QLabel(tr("[%1x%2, %3 Mpx]").arg(width).arg(height).arg(mpxs, 0, 'f', 2));
    QLabel* sizeLabel = new QLabel(tr("[%1x%2]").arg(width).arg(height));
    QToolButton* exportButton = new QToolButton();
    exportButton->setText(tr("Export..."));

    QHBoxLayout* bottomLayout = new QHBoxLayout();
    bottomLayout->addWidget(viewLabel);
    bottomLayout->addWidget(originalRadio);
    bottomLayout->addWidget(processRadio);
    bottomLayout->addWidget(new QLabel(tr("Zoom:")));
    bottomLayout->addWidget(zoomLabel);
    bottomLayout->addWidget(fullButton);
    bottomLayout->addWidget(fitButton);
    bottomLayout->addStretch();
    bottomLayout->addWidget(exportButton);
    //bottomLayout->addStretch();
    bottomLayout->addWidget(sizeLabel);

    QVBoxLayout* mainLayout = new QVBoxLayout();
    if (!title.isEmpty())
    {
        QLabel* titleLabel = new QLabel(title);
        Utility::changeFont(titleLabel, true);
        titleLabel->setAlignment(Qt::AlignCenter);
        mainLayout->addWidget(titleLabel);
    }
    mainLayout->addWidget(imageView);
    mainLayout->addLayout(bottomLayout);
    setLayout(mainLayout);

    connect(originalRadio, &QRadioButton::toggled,
            this,          &SingleViewer::toggleMode);
    connect(fitButton,  &QToolButton::clicked,
            imageView,  &DynamicView::fitZoom);
    connect(fullButton, &QToolButton::clicked,
            imageView,  &DynamicView::fullZoom);
    connect(exportButton, &QToolButton::clicked,
            this,         &SingleViewer::exportImage);
    connect(imageView,  &DynamicView::viewChanged,
            this,       &SingleViewer::forwardViewChanged);

    viewLabel->setVisible(!processed.empty());
    originalRadio->setVisible(!processed.empty());
    processRadio->setVisible(!processed.empty());
    exportButton->setVisible(!processed.empty());
    if (!processed.empty())
    {
        originalRadio->setChecked(false);
        processRadio->setChecked(true);
        toggleMode(false);
    }
}

void SingleViewer::updateProcessed(const cv::Mat &processed)
{
    if (this->processed.empty())
        return;
    this->processed = processed;
    toggleMode(originalRadio->isChecked());
}

void SingleViewer::updateOriginal(const cv::Mat &original)
{
    this->original = original;
    toggleMode(true);
}

void SingleViewer::changeView(QRect, double zoomFactor, int horizScroll, int vertScroll)
{
    imageView->changeView(QRect(), zoomFactor, horizScroll, vertScroll);
}

void SingleViewer::forwardViewChanged(QRect sceneRect, double zoomFactor, int horizScroll, int vertScroll)
{
    QString text = QString::number(zoomFactor * 100, 'f', 2) + "%";
    /*
    if (zoom_level < 100)
        text += tr(" [missing pixels]");
    else if (zoom_level > 100)
        text += tr(" [interpolated pixels]");
    else
        text += tr(" [no interpolation]");
    */
    zoomLabel->setText(text);
    Utility::changeFont(zoomLabel, zoomFactor == 1);
    emit viewChanged(sceneRect, zoomFactor, horizScroll, vertScroll);
}

QRect SingleViewer::getSceneRect()
{
    return imageView->getSceneRect();
}

void SingleViewer::keyPressEvent(QKeyEvent* event)
{
    if (event->key() == Qt::Key_Space)
    {
        if (originalRadio->isChecked())
            processRadio->setChecked(true);
        else
            originalRadio->setChecked(true);
    }
    else
        QWidget::keyPressEvent(event);
}

void SingleViewer::toggleMode(bool origToggled)
{
    if (origToggled)
        imageItem->setPixmap(QPixmap::fromImage(QtOcv::mat2Image(original)));
    else if (!processed.empty())
        imageItem->setPixmap(QPixmap::fromImage(QtOcv::mat2Image(processed)));
}

void SingleViewer::exportImage()
{
    QString fileName = QFileDialog::getSaveFileName(this, QApplication::applicationName(),
                                                    QString(), tr("PNG Images (*.png)"));
    if (fileName.isEmpty())
        return;
    if (QFileInfo(fileName).suffix() != "png")
        fileName.append(".png");
    cv::Mat exported = originalRadio->isChecked() ? original : processed;
    cv::imwrite(fileName.toStdString(), exported);
}
