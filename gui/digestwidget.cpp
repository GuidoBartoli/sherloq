#include "digestwidget.h"
#include "utility.h"

DigestWidget::DigestWidget(const QString &fileName, QWidget *parent) : ToolWidget(parent)
{
    //const QString TIME_FORMAT = QString("yyyy-MM-dd HH:mm:ss.z t");
    const int ROWS = 8;
    const int COLUMNS = 2;

    QTableWidget* tableWidget = new QTableWidget();
    tableWidget->setColumnCount(COLUMNS);
    tableWidget->setRowCount(ROWS);
    QStringList headers;
    headers << tr("Property") << tr("Value");
    tableWidget->setHorizontalHeaderLabels(headers);

    QFileInfo file_info(fileName);
    tableWidget->setItem(0, 0, new QTableWidgetItem(tr("File name")));
    tableWidget->setItem(0, 1, new QTableWidgetItem(file_info.fileName()));
    QLocale locale;
    tableWidget->setItem(1, 0, new QTableWidgetItem(tr("Physical size")));
    tableWidget->setItem(1, 1, new QTableWidgetItem(locale.toString(file_info.size()) + " bytes"));
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly))
        return;

    QByteArray bytes = file.readAll();
    QString md5(QCryptographicHash::hash(bytes, QCryptographicHash::Md5).toHex());
    tableWidget->setItem(2, 0, new QTableWidgetItem(tr("MD5")));
    tableWidget->setItem(2, 1, new QTableWidgetItem(md5));
    QString sha1(QCryptographicHash::hash(bytes, QCryptographicHash::Sha1).toHex());
    tableWidget->setItem(3, 0, new QTableWidgetItem(tr("SHA-1")));
    tableWidget->setItem(3, 1, new QTableWidgetItem(sha1));
    QString sha256(QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex());
    tableWidget->setItem(4, 0, new QTableWidgetItem(tr("SHA-256")));
    tableWidget->setItem(4, 1, new QTableWidgetItem(sha256));

    QString byteString;
    QString bitString;
    computePHash(fileName, byteString, bitString);
    tableWidget->setItem(5, 0, new QTableWidgetItem(tr("pHash (HEX)")));
    tableWidget->setItem(5, 1, new QTableWidgetItem(byteString));
    tableWidget->setItem(6, 0, new QTableWidgetItem(tr("pHash (binary)")));
    tableWidget->setItem(6, 1, new QTableWidgetItem(bitString));

    QString device = nameBallistics(file_info.fileName());
    tableWidget->setItem(7, 0, new QTableWidgetItem(tr("Name ballistics")));
    tableWidget->setItem(7, 1, new QTableWidgetItem(device.isEmpty() ? tr("Unknown or manually renamed") : device));

    for (int i = 0; i < ROWS; ++i)
    {
        QFont font = tableWidget->item(i, 0)->font();
        font.setBold(true);
        tableWidget->item(i, 0)->setFont(font);
    }
    tableWidget->resizeColumnsToContents();
    tableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);

    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addWidget(tableWidget);
    setLayout(vertLayout);
    setMinimumSize(750, 300);

    connect(tableWidget, &QTableWidget::itemDoubleClicked,
            this,        &DigestWidget::copyCell);
}

void DigestWidget::copyCell(QTableWidgetItem* item)
{
    QApplication::clipboard()->setText(item->text());
    emit messageToShow(tr("[GENERAL::Digest] Cell contents copied to clipboard"));
}

void DigestWidget::computePHash(const QString &fileName, QString &bytes, QString &bits)
{
    cv::Mat image = cv::imread(fileName.toStdString(), cv::IMREAD_GRAYSCALE);
    cv::Mat filtered;
    cv::blur(image, filtered, cv::Size(7, 7));
    cv::Mat block;
    cv::resize(filtered, block, cv::Size(32, 32), cv::INTER_AREA);
    block.convertTo(block, CV_32F, 1, -128);
    cv::Mat coeffs;
    cv::dct(block, coeffs);

    cv::Mat low(coeffs, cv::Rect(0, 0, 8, 8));
    cv::Mat row = low.clone().reshape(0, 1);
    std::vector<double> vec;
    row.copyTo(vec);
    std::nth_element(vec.begin(), vec.begin() + vec.size() / 2, vec.end());
    double median = vec.at(vec.size() / 2);

    quint64 one = 0x1000000000000000;
    quint64 hash = 0x0000000000000000;
    for (int k = 0; k < row.cols; ++k)
    {
        int i = Utility::ZIG_ZAG.at<int>(k, 0);
        int j = Utility::ZIG_ZAG.at<int>(k, 1);
        float coeff = coeffs.at<float>(i, j);
        if (coeff > median)
            hash |= one;
        one >>= 1;
    }

    /*
    quint64 one = 0x0000000000000001;
    quint64 hash = 0x0000000000000000;
    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 8; ++j)
        {
            if (low.at<float>(i, j) > median)
                hash |= one;
            one = one << 1;
        }
    }
    */

    bytes = QString(QByteArray::number(hash, 16));
    bits = QString(QByteArray::number(hash, 2));
}

QString DigestWidget::nameBallistics(const QString &fileName)
{
    QList< QPair<QString,QString> > patterns;
    patterns.append(qMakePair(tr("^DSCN[0-9]{4}\\.JPG$"), tr("Nikon Coolpix")));
    patterns.append(qMakePair(tr("^DSC_[0-9]{4}\\.JPG$"), tr("Nikon DSLR")));
    patterns.append(qMakePair(tr("^IMG_[0-9]{4}\\.JPG$"), tr("Canon camera")));
    patterns.append(qMakePair(tr("^IMAG[0-9]{4}\\.JPG$"), tr("HTC One")));
    patterns.append(qMakePair(tr("^FUJI[0-9]{4}\\.JPG$"), tr("Fujifilm camera")));

    QRegularExpression::PatternOptions opts = QRegularExpression::CaseInsensitiveOption;
    for (int i = 0; i < patterns.length(); ++i)
        if (QRegularExpression(patterns.at(i).first, opts).match(fileName).hasMatch())
            return patterns.at(i).second;
    return QString();
}

/*
double DigestWidget::medianMat(cv::Mat Input)
{
    Input = Input.reshape(0,1); // spread Input Mat to single row
    std::vector<double> vecFromMat;
    Input.copyTo(vecFromMat); // Copy Input Mat to vector vecFromMat
    vecFromMat.erase(vecFromMat.begin());
    //qDebug() << vecFromMat.size();
    std::nth_element(vecFromMat.begin(), vecFromMat.begin() + vecFromMat.size() / 2, vecFromMat.end());
    return vecFromMat[vecFromMat.size() / 2];
}
*/

