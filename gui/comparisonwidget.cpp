#include "comparisonwidget.h"
#include "utility.h"

ComparisonWidget::ComparisonWidget(const cv::Mat &evidence, QWidget *parent)
    : ToolWidget(parent),
      evidence(evidence)
{
    filenameLabel = new QLabel();
    Utility::changeFont(filenameLabel, true);
    indicesLabel = new QLabel();

    QPushButton* loadButton = new QPushButton(tr("Load reference..."));
    connect(loadButton, &QPushButton::clicked,
            this,       &ComparisonWidget::loadReference);
    originalRadio = new QRadioButton(tr("Original"));
    differenceRadio = new QRadioButton(tr("Difference"));
    indexMapRadio = new QRadioButton(tr("Index map"));
    originalRadio->setChecked(true);
    connect(originalRadio, &QRadioButton::toggled,
            this,          &ComparisonWidget::updateReference);
    connect(differenceRadio, &QRadioButton::toggled,
            this,            &ComparisonWidget::updateReference);
    connect(indexMapRadio, &QRadioButton::toggled,
            this,          &ComparisonWidget::updateReference);

    QHBoxLayout* horizLayout = new QHBoxLayout();
    horizLayout->addWidget(loadButton);
    horizLayout->addWidget(filenameLabel);
    horizLayout->addWidget(indicesLabel);
    horizLayout->addStretch();
    horizLayout->addWidget(originalRadio);
    horizLayout->addWidget(differenceRadio);
    horizLayout->addWidget(indexMapRadio);

    evidenceViewer = new SingleViewer(evidence, cv::Mat(), tr("Evidence"));
    cv::Mat gray(evidence.rows, evidence.cols, CV_8UC3, cv::Scalar(127, 127, 127));
    referenceViewer = new SingleViewer(gray, cv::Mat(), tr("Reference"));
    connect(evidenceViewer, &SingleViewer::viewChanged,
            referenceViewer, &SingleViewer::changeView);
    connect(referenceViewer, &SingleViewer::viewChanged,
            evidenceViewer,  &SingleViewer::changeView);
    QHBoxLayout* viewLayout = new QHBoxLayout();
    viewLayout->addWidget(evidenceViewer);
    viewLayout->addWidget(referenceViewer);

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addLayout(horizLayout);
    vertLayout->addLayout(viewLayout);
    setLayout(vertLayout);
}

void ComparisonWidget::loadReference()
{
    QSettings settings;
    QString fileName = QFileDialog::getOpenFileName(
                this, tr("Load image"), settings.value("last_folder_ref").toString(),
                tr("Supported images (*.jpg *.jpeg *.png *.tif *.tiff)"));
    if (fileName.isEmpty())
        return;
    cv::Mat image = cv::imread(fileName.toStdString(), CV_LOAD_IMAGE_COLOR);
    if (image.empty())
    {
        QMessageBox::critical(this, QApplication::applicationName(), tr("Unable to load reference image!"));
        return;
    }
    if (image.channels() > 3)
    {
        QMessageBox::warning(this, QApplication::applicationName(), tr("Embedded alpha channel discarded"));
        image.convertTo(image, CV_BGRA2BGR);
    }
    if (image.size != evidence.size)
    {
        QMessageBox::critical(this, QApplication::applicationName(), tr("Evidence and Reference images must have the same size"));
        return;
    }
    reference = image;
    difference = cv::abs(evidence - reference);
    cv::cvtColor(difference, difference, CV_BGR2GRAY);
    cv::normalize(difference, difference, 0, 255, cv::NORM_MINMAX);
    cv::cvtColor(difference, difference, CV_GRAY2BGR);
    referenceViewer->updateOriginal(reference);

    double psnrValue = psnr(evidence, reference);
    double ssimValue = ssim2(evidence, reference);
    QString psnrString = psnrValue == -1 ? tr("+inf") : QString::number(psnrValue, 'f', 4);
    QString ssimString = QString::number(ssimValue, 'f', 4);
    indicesLabel->setText(tr("(PSNR = %1 dB, SSIM = %2)").arg(psnrString).arg(ssimString));
    originalRadio->setChecked(true);

    filenameLabel->setText(QFileInfo(fileName).fileName());
}

void ComparisonWidget::updateReference()
{
    if (originalRadio->isChecked())
        referenceViewer->updateOriginal(reference);
    else if (differenceRadio->isChecked())
        referenceViewer->updateOriginal(difference);
    else if (indexMapRadio->isChecked())
        referenceViewer->updateOriginal(indexMap);
}

// ------------------------------------------------------------------------------------------------

double ComparisonWidget::sigma(const cv::Mat & m, int i, int j, int block_size)
{
    double sd = 0;
    cv::Mat m_tmp = m(cv::Range(i, i + block_size), cv::Range(j, j + block_size));
    cv::Mat m_squared(block_size, block_size, CV_64F);
    cv::multiply(m_tmp, m_tmp, m_squared);

    // E(x)
    double avg = cv::mean(m_tmp)[0];
    // E(x²)
    double avg_2 = cv::mean(m_squared)[0];
    sd = sqrt(avg_2 - avg * avg);
    return sd;
}

// Covariance
double ComparisonWidget::cov(const cv::Mat &m1, const cv::Mat &m2, int i, int j, int block_size)
{
    cv::Mat m3 = cv::Mat::zeros(block_size, block_size, m1.depth());
    cv::Mat m1_tmp = m1(cv::Range(i, i + block_size), cv::Range(j, j + block_size));
    cv::Mat m2_tmp = m2(cv::Range(i, i + block_size), cv::Range(j, j + block_size));
    cv::multiply(m1_tmp, m2_tmp, m3);

    double avg_ro 	= mean(m3)[0]; // E(XY)
    double avg_r 	= mean(m1_tmp)[0]; // E(X)
    double avg_o 	= mean(m2_tmp)[0]; // E(Y)
    double sd_ro = avg_ro - avg_o * avg_r; // E(XY) - E(X)E(Y)
    return sd_ro;
}

// Mean squared error
double ComparisonWidget::eqm(const cv::Mat &img1, const cv::Mat &img2)
{
    int i, j;
    double eqm = 0;
    int height = img1.rows;
    int width = img1.cols;

    for (i = 0; i < height; i++)
        for (j = 0; j < width; j++)
            eqm += (img1.at<double>(i, j) - img2.at<double>(i, j)) * (img1.at<double>(i, j) - img2.at<double>(i, j));

    eqm /= height * width;

    return eqm;
}

double ComparisonWidget::psnr(const cv::Mat &img1, const cv::Mat &img2)
{
    cv::Mat img1gray, img2gray;
    cv::cvtColor(img1, img1gray, CV_BGR2GRAY);
    cv::cvtColor(img2, img2gray, CV_BGR2GRAY);
    img1gray.convertTo(img1gray, CV_32F);
    img2gray.convertTo(img2gray, CV_32F);
    cv::Mat diff;
    cv::pow(img1gray - img2gray, 2, diff);
    double den = sqrt(cv::mean(diff)[0]);
    if (den == 0)
        return -1;
    return (20 * log10((255 * 255) / den));
}

double ComparisonWidget::ssim(const cv::Mat &img_src, const cv::Mat &img_compressed, int block_size)
{
    const float C1 = (0.01 * 255 * 0.01  * 255);
    const float C2 = (0.03 * 255 * 0.03  * 255);
    double ssim = 0;

    //TODO: Aggiungere padImage per considerare tutta l'immagine

    int nbBlockPerHeight = img_src.rows / block_size;
    int nbBlockPerWidth  = img_src.cols / block_size;

    cv::Mat img1gray, img2gray;
    cv::cvtColor(img_src, img1gray, CV_BGR2GRAY);
    cv::cvtColor(img_compressed, img2gray, CV_BGR2GRAY);
    img1gray.convertTo(img1gray, CV_32F);
    img2gray.convertTo(img2gray, CV_32F);

    for (int k = 0; k < nbBlockPerHeight; k++)
    {
        for (int l = 0; l < nbBlockPerWidth; l++)
        {
            int m = k * block_size;
            int n = l * block_size;

            double avg_o 	= cv::mean(img1gray(cv::Range(k, k + block_size), cv::Range(l, l + block_size)))[0];
            double avg_r 	= cv::mean(img2gray(cv::Range(k, k + block_size), cv::Range(l, l + block_size)))[0];
            double sigma_o 	= sigma(img1gray, m, n, block_size);
            double sigma_r 	= sigma(img2gray, m, n, block_size);
            double sigma_ro	= cov(img1gray, img2gray, m, n, block_size);

            ssim += ((2 * avg_o * avg_r + C1) * (2 * sigma_ro + C2)) / ((avg_o * avg_o + avg_r * avg_r + C1) * (sigma_o * sigma_o + sigma_r * sigma_r + C2));

        }
    }
    ssim /= nbBlockPerHeight * nbBlockPerWidth;

    return ssim;
}

double ComparisonWidget::ssim2(const cv::Mat &img1orig, const cv::Mat &img2orig)
{
    const double C1 = 6.5025;
    const double C2 = 58.5225;
    const int    K  = 11;
    const double S  = 1.5;

    cv::Mat img1, img2, img1img2;
    cv::cvtColor(img1orig, img1, CV_BGR2GRAY);
    cv::cvtColor(img2orig, img2, CV_BGR2GRAY);
    img1.convertTo(img1, CV_32F);
    img2.convertTo(img2, CV_32F);
    cv::Mat img1sqr, img2sqr;
    cv::pow(img1, 2, img1sqr);
    cv::pow(img2, 2, img2sqr);
    cv::multiply(img1, img2, img1img2);

    cv::Mat mu1, mu2;
    cv::GaussianBlur(img1, mu1, cv::Size(K, K), S);
    cv::GaussianBlur(img2, mu2, cv::Size(K, K), S);
    cv::Mat mu1sqr, mu2sqr, mu1mu2;
    cv::pow(mu1, 2, mu1sqr);
    cv::pow(mu2, 2, mu2sqr);
    cv::multiply(mu1, mu2, mu1mu2);

    cv::Mat sig1sqr, sig2sqr, sig1sig2;
    cv::GaussianBlur(img1sqr, sig1sqr, cv::Size(K, K), S);
    cv::addWeighted(sig1sqr, 1, mu1sqr, -1, 0, sig1sqr);
    cv::GaussianBlur(img2sqr, sig2sqr, cv::Size(K, K), S);
    cv::addWeighted(sig2sqr, 1, mu2sqr, -1, 0, sig2sqr);
    cv::GaussianBlur(img1img2, sig1sig2, cv::Size(K, K), S);
    cv::addWeighted(sig1sig2, 1, mu1mu2, -1, 0, sig1sig2);

    cv::Mat num = (2*mu1mu2 + C1).mul(2*sig1sig2 + C2);
    cv::Mat den = (mu1sqr + mu2sqr + C1).mul(sig1sqr + sig2sqr + C2);
    indexMap = num / den;
    double ssim = cv::mean(indexMap)[0];

    cv::normalize(indexMap, indexMap, 0, 255, cv::NORM_MINMAX);
    indexMap.convertTo(indexMap, CV_8U);
    cv::cvtColor(indexMap, indexMap, CV_GRAY2BGR);

    return ssim;
}
