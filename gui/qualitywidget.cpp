#include "qualitywidget.h"
#include "ExifTool.h"
#include "utility.h"

QualityWidget::QualityWidget(const QString &fileName, QWidget* parent)
    : ToolWidget(parent)
{
    const QString TEMP_FILE = QDir::tempPath() + "/quality.jpg";

    const uchar JPEG_MRK = 0xFF;
    const uchar JPEG_SOI = 0xD8;
    const uchar JPEG_DQT = 0xDB;
    const uchar JPEG_PAD = 0x00;

    const int MAX_TABLES = 2;
    const int LEN_OFFSET = 2;
    const int START_QUAL = 0;
    const int END_QUAL   = 100;
    const int QUAL_STEP  = 1;
    const int NUM_STEPS  = (END_QUAL - START_QUAL) / QUAL_STEP + 1;
    const int LUMA_IDX   = 0;
    const int CHROMA_IDX = 1;

    QFile(fileName).copy(TEMP_FILE);
    ExifTool* exif = new ExifTool();
    exif->SetNewValue("all", NULL);
    exif->WriteInfo(TEMP_FILE.toStdString().c_str());
    int result = exif->Complete();
    delete exif;
    if (result <= 0)
    {
        QFile(TEMP_FILE).remove();
        showError(0);
        return;
    }

    FILE* file = fopen(TEMP_FILE.toStdString().c_str(), "rb");
    if (fgetc(file) != JPEG_MRK || fgetc(file) != JPEG_SOI)
    {
        showError(1);
        fclose(file);
        QFile(TEMP_FILE).remove();
        return;
    }

    bool found = false;
    cv::Mat luma(Utility::DCT_SIZE, Utility::DCT_SIZE, CV_8U);
    cv::Mat chroma(Utility::DCT_SIZE, Utility::DCT_SIZE, CV_8U);
    while (!feof(file))
    {
        if (!findNextMarker(file, JPEG_MRK, JPEG_DQT, JPEG_PAD))
            break;
        int length = fgetc(file) - LEN_OFFSET;
        if ((length % (Utility::TABLE_SIZE + 1)) != 0 || length <= 0)
            continue;
        while (length > 0)
        {
            int type = fgetc(file);
            --length;
            int index = type & 0x0f;
            if (index >= MAX_TABLES)
                break;
            for (int k = 0; k < Utility::TABLE_SIZE; ++k)
            {
                uchar c = fgetc(file);
                --length;
                if (feof(file))
                    break;
                int i = Utility::ZIG_ZAG.at<int>(k, 0);
                int j = Utility::ZIG_ZAG.at<int>(k, 1);
                if (index == LUMA_IDX)
                    luma.at<uchar>(i, j) = c;
                else if (index == CHROMA_IDX)
                    chroma.at<uchar>(i, j) = c;
            }
            if (luma.at<uchar>(0, 0) != 0)
                found = true;
        }
    }
    fclose(file);
    QFile(TEMP_FILE).remove();
    if (!found)
    {
        showError(2);
        return;
    }

    float luma_avg = mean(luma)[0];
    luma_avg *= Utility::TABLE_SIZE;
    luma_avg -= luma.at<uchar>(0, 0);
    luma_avg /= Utility::TABLE_SIZE - 1;
    float luma_q = (1 - luma_avg / 255) * 100;
    float chroma_avg = mean(chroma)[0];
    chroma_avg *= Utility::TABLE_SIZE;
    chroma_avg -= chroma.at<uchar>(0, 0);
    chroma_avg /= Utility::TABLE_SIZE - 1;
    float chroma_q = (1 - chroma_avg / 255) * 100;
    cv::Mat search(1, NUM_STEPS, CV_32F, cv::Scalar(0));
    int count = 0;
    for (int q = START_QUAL; q <= END_QUAL; q += QUAL_STEP)
    {
        cv::Mat luma0;
        cv::Mat chroma0;
        getTables(q, luma0, chroma0);
        float diff = mean(abs(luma - luma0))[0] + mean(abs(chroma - chroma0))[0];
        search.at<float>(0, count++) = diff;
    }
    double min_val;
    cv::Point min_loc;
    minMaxLoc(search, &min_val, 0, &min_loc);
    int quality = START_QUAL + QUAL_STEP * min_loc.x;


    //fixedFont.setPointSize(fixedFont.pointSize() - 1);

    QPlainTextEdit* lumaEdit = new QPlainTextEdit();
    //lumaEdit->setFont(Utility::FIXED_FONT);
    lumaEdit->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    lumaEdit->setReadOnly(true);
    lumaEdit->appendPlainText(tr("Luminance QT (level = %1%)\n").arg(luma_q, 0, 'f', 4));
    for (int i = 0; i < Utility::DCT_SIZE; ++i)
    {
        QString line;
        for (int j = 0; j < Utility::DCT_SIZE; ++j)
        {
            uchar lumaValue = luma.at<uchar>(i, j);
            QString text = QString::number(lumaValue);
            if (lumaValue < 10)
                text = "  " + text;
            else if (lumaValue < 100)
                text = " " + text;
            line.append(text + " ");
        }
        lumaEdit->appendPlainText(line);
    }
    lumaEdit->setMinimumSize(320, 220);

    QPlainTextEdit* chromaEdit = new QPlainTextEdit();
    chromaEdit->setFont(lumaEdit->font());
    chromaEdit->setReadOnly(true);
    chromaEdit->appendPlainText(tr("Chrominance QT (level = %1%)\n").arg(chroma_q, 0, 'f', 4));
    for (int i = 0; i < Utility::DCT_SIZE; ++i)
    {
        QString line;
        for (int j = 0; j < Utility::DCT_SIZE; ++j)
        {
            uchar chromaValue = chroma.at<uchar>(i, j);
            QString text = QString::number(chromaValue);
            if (chromaValue < 10)
                text = "  " + text;
            else if (chromaValue < 100)
                text = " " + text;
            line.append(text + " ");
        }
        chromaEdit->appendPlainText(line);
    }
    chromaEdit->setMinimumSize(lumaEdit->minimumSize());

    QLabel* qualityLabel = new QLabel(tr("Estimated JPEG quality (last save) = %1%").arg(quality));
    QLabel* deviationLabel = new QLabel(tr("(deviation from standard QTs = %2)").arg(min_val, 0, 'f', 2));
    deviationLabel->setAlignment(Qt::AlignRight);
    Utility::changeFont(qualityLabel, true);
    Utility::changeFont(deviationLabel, false, true);

    QGridLayout* gridLayout = new QGridLayout();
    gridLayout->addWidget(lumaEdit, 0, 0);
    gridLayout->addWidget(chromaEdit, 0, 1);
    gridLayout->addWidget(qualityLabel, 1, 0);
    gridLayout->addWidget(deviationLabel, 1, 1);
    setLayout(gridLayout);
}

void QualityWidget::getTables(int quality, cv::Mat &luma_qt, cv::Mat &chroma_qt, bool baseline)
{
    const cv::Mat LUMA_QT = (cv::Mat_<double>(Utility::DCT_SIZE, Utility::DCT_SIZE) <<
        16,  11,  10,  16,  24,  40,  51,  61,
        12,  12,  14,  19,  26,  58,  60,  55,
        14,  13,  16,  24,  40,  57,  69,  56,
        14,  17,  22,  29,  51,  87,  80,  62,
        18,  22,  37,  56,  68, 109, 103,  77,
        24,  35,  55,  64,  81, 104, 113,  92,
        49,  64,  78,  87, 103, 121, 120, 101,
        72,  92,  95,  98, 112, 100, 103,  99);
    const cv::Mat CHROMA_QT = (cv::Mat_<double>(Utility::DCT_SIZE, Utility::DCT_SIZE) <<
        17,  18,  24,  47,  99,  99,  99,  99,
        18,  21,  26,  66,  99,  99,  99,  99,
        24,  26,  56,  99,  99,  99,  99,  99,
        47,  66,  99,  99,  99,  99,  99,  99,
        99,  99,  99,  99,  99,  99,  99,  99,
        99,  99,  99,  99,  99,  99,  99,  99,
        99,  99,  99,  99,  99,  99,  99,  99,
        99,  99,  99,  99,  99,  99,  99,  99);

    if (quality <= 0)
        quality = 1;
    else if (quality > 100)
        quality = 100;
    if (quality < 50)
        quality = 5000 / quality;
    else
        quality = 200 - quality * 2;

    luma_qt = (LUMA_QT * quality + 50) / 100;
    luma_qt.convertTo(luma_qt, baseline ? CV_8U : CV_16U);
    luma_qt.setTo(1, luma_qt <= 0);
    chroma_qt = (CHROMA_QT * quality + 50) / 100;
    chroma_qt.convertTo(chroma_qt, baseline ? CV_8U : CV_16U);
    chroma_qt.setTo(1, chroma_qt <= 0);
}

void QualityWidget::showError(int type)
{
    QLabel* errorLabel = new QLabel();
    if (type == 0)
        errorLabel->setText(tr("Error while reading file!"));
    else if (type == 1)
        errorLabel->setText(tr("File is not a JPEG image!"));
    else if (type == 2)
        errorLabel->setText(tr("Unable to find JPEG tables!"));
    else
        errorLabel->setText(tr("Unknown error!"));
    errorLabel->setAlignment(Qt::AlignCenter);
    QFont font = errorLabel->font();
    font.setItalic(true);
    errorLabel->setFont(font);
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addWidget(errorLabel);
    setLayout(vertLayout);
}

bool QualityWidget::findNextMarker(FILE* file, uchar byte1, uchar byte2, uchar byte3)
{
    while (!feof(file))
        if (fgetc(file) == byte1 && fgetc(file) == byte2 && fgetc(file) == byte3)
            return true;
    return false;
}
