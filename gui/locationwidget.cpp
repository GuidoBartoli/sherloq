#include "locationwidget.h"
#include <QtWebEngineWidgets>

LocationWidget::LocationWidget(const QString &fileName, QWidget* parent) : ToolWidget(parent)
{
    const int   GEO_PRECISION = 6;
    const QChar GEO_SEPARATOR = QChar(32);

    ExifTool* exif = new ExifTool();
    TagInfo* info = exif->ImageInfo(fileName.toStdString().c_str(), "-c\n\"%+.6f\"\n-GPSPosition");
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
    QString position = QString::fromLocal8Bit(info->next->num, info->next->numLen);
    QStringList latlong = position.split(GEO_SEPARATOR);
    double lat = latlong.at(0).toDouble();
    double lon = latlong.at(1).toDouble();
    QUrl url(tr("http://maps.google.com/maps?z=12&t=m&q=loc:%1+%2").arg(lat, 0, 'f', GEO_PRECISION).arg(lon, 0, 'f', GEO_PRECISION));
    //QUrl url("http://maps.google.co.uk/maps/place/42.3525,-71.0912/@42.3525,-71.0912,12z/data=!3m1!1e3");
    //qDebug() << url;

    QWebEngineView* webView = new QWebEngineView();
    webView->load(url);
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addWidget(webView);
    setLayout(vertLayout);

//    http://maps.google.com/maps?z=12&t=m&q=loc:38.9419+-78.3020
//    - z is the zoom level (1-20)
//    - t is the map type ("m" map, "k" satellite, "h" hybrid, "p" terrain, "e" GoogleEarth)
//    - q is the search query, if it is prefixed by loc: then google assumes it is a lat lon separated by a +

    // http://maps.google.co.uk/maps/place/52.03877,-2.3416/@52.03877,-2.3416,15z/data=!3m1!1e3
    //
}

void LocationWidget::showError(int type, ExifTool *exif, TagInfo *info)
{
    QLabel* errorLabel = new QLabel();
    if (type == 0)
        errorLabel->setText(tr("Error while extracting information from file!"));
    else if (type == 1)
        errorLabel->setText(tr("File does not contain any geo-location data"));
    else if (type == 2)
        errorLabel->setText(tr("Error while reading geo-location data"));
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
