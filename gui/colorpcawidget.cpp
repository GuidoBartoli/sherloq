#include "colorpcawidget.h"

ColorPcaWidget::ColorPcaWidget(const cv::Mat &input, QWidget* parent) : ToolWidget(parent)
{
    const int DIMS = 3;

    QProgressDialog progress(tr("Computing PCA..."), QString(), 0, 1 + DIMS*5, parent);
    progress.setWindowModality(Qt::WindowModal);

    cv::Mat rgbData(input.rows * input.cols, DIMS, CV_32F);
    cv::MatConstIterator_<cv::Vec3b> inputIt = input.begin<cv::Vec3b>();
    cv::MatConstIterator_<cv::Vec3b> inputEnd = input.end<cv::Vec3b>();
    int i = 0;
    while (inputIt != inputEnd)
    {
        for (int j = 0; j < DIMS; ++j)
            rgbData.at<float>(i, j) = (*inputIt)[j];
        ++inputIt;
        ++i;
    }
    cv::PCA pca(rgbData, cv::noArray(), CV_PCA_DATA_AS_ROW, DIMS);
    progress.setValue(1);

    cv::Point3f q(pca.mean.at<float>(0), pca.mean.at<float>(1), pca.mean.at<float>(2));
    for (int d = 0; d < DIMS; ++d)
    {
        cv::Point3f v = cv::Point3f(pca.eigenvectors.at<float>(d, 0),
                                    pca.eigenvectors.at<float>(d, 1),
                                    pca.eigenvectors.at<float>(d, 2));
        cv::Point3f r = q + v * pca.eigenvalues.at<float>(d);
        cv::Mat distance(input.rows, input.cols, CV_32F);
        cv::Mat cross(input.rows, input.cols, CV_32FC3);
        cv::MatIterator_<float> distance_it = distance.begin<float>();
        cv::MatIterator_<cv::Vec3f> cross_it = cross.begin<cv::Vec3f>();
        inputIt = input.begin<cv::Vec3b>();
        while (inputIt != inputEnd)
        {
            cv::Point3f p((*inputIt)[0], (*inputIt)[1], (*inputIt)[2]);
            *distance_it = cv::norm((p - q).cross(p - r)) / cv::norm(r - q);
            *cross_it = p.cross(v);
            ++inputIt;
            ++distance_it;
            ++cross_it;
        }
        progress.setValue(progress.value() + 1);

        std::vector<cv::Mat> channels;
        cv::split(cross, channels);
        for (int c = 0; c < 3; ++c)
        {
            cv::Mat channel;
            cv::normalize(channels.at(c), channel, 0, 255, cv::NORM_MINMAX);
            channel.convertTo(channel, CV_8U);
            channels[c] = channel;
        }
        progress.setValue(progress.value() + 1);

        cv::normalize(distance, distance, 0, 255, cv::NORM_MINMAX);
        distance.convertTo(distance, CV_8U);
        components.push_back(distance);
        progress.setValue(progress.value() + 1);

        cv::Mat distance_eq;
        cv::equalizeHist(distance, distance_eq);
        components.push_back(distance_eq);
        progress.setValue(progress.value() + 1);

        cv::Mat cross_img;
        cv::merge(channels, cross_img);
        components.push_back(cross_img);
        progress.setValue(progress.value() + 1);
    }

    componentCombo = new QComboBox();
    componentCombo->addItem(tr("First"));
    componentCombo->addItem(tr("Second"));
    componentCombo->addItem(tr("Third"));
    distanceRadio = new QRadioButton(tr("Vector Distance"));
    distanceEqRadio = new QRadioButton(tr("Distance Equalized"));
    crossRadio = new QRadioButton(tr("Cross-Correlation"));
    distanceRadio->setChecked(true);
    lastRadio = distanceRadio;

    connect(componentCombo, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
            this,           &ColorPcaWidget::updateView);
    connect(distanceRadio, &QRadioButton::toggled,
            this,          &ColorPcaWidget::updateView);
    connect(distanceEqRadio, &QRadioButton::toggled,
            this,            &ColorPcaWidget::updateView);
    connect(crossRadio, &QRadioButton::toggled,
            this,       &ColorPcaWidget::updateView);

    QHBoxLayout* horizLayout = new QHBoxLayout();
    horizLayout->addWidget(new QLabel(tr("Component:")));
    horizLayout->addWidget(componentCombo);
    horizLayout->addWidget(distanceRadio);
    horizLayout->addWidget(distanceEqRadio);
    horizLayout->addWidget(crossRadio);
    horizLayout->addStretch();

    pcaViewer = new SingleViewer(input, components.at(0));

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(horizLayout);
    vertLayout->addWidget(pcaViewer);
    setLayout(vertLayout);
}

void ColorPcaWidget::updateView()
{
    int index = componentCombo->currentIndex();
    if (distanceRadio->isChecked())
    {
        pcaViewer->updateProcessed(components.at(index * 3));
        lastRadio = distanceRadio;
    }
    else if (distanceEqRadio->isChecked())
    {
        pcaViewer->updateProcessed(components.at(index * 3 + 1));
        lastRadio = distanceEqRadio;
    }
    else if (crossRadio->isChecked())
    {
        pcaViewer->updateProcessed(components.at(index * 3 + 2));
        lastRadio = crossRadio;
    }
    else
        lastRadio->setChecked(true);
}
