#include "structurewidget.h"
#include <QtWebEngineWidgets>

StructureWidget::StructureWidget(const QString &fileName, QWidget* parent) : ToolWidget(parent)
{
    QProcess process;
    process.start("exiftool -htmldump0 " + fileName);
    if (!process.waitForFinished())
    {
        emit messageToShow(tr("[FILE::Structure] Unable to extract data!"));
        return;
    }
    QByteArray bytes = process.readAllStandardOutput();
    QFile file(TEMP_FILE);
    if (!file.open(QIODevice::WriteOnly))
    {
        emit messageToShow(tr("[FILE::Structure] Unable to dump file!"));
        return;
    }
    file.write(bytes);
    file.close();

    QWebEngineView* webView = new QWebEngineView();
    webView->load(QUrl("file://" + TEMP_FILE));
    QVBoxLayout* vertLayout = new QVBoxLayout();
    vertLayout->addWidget(webView);
    setLayout(vertLayout);
}

StructureWidget::~StructureWidget()
{
    QFile(TEMP_FILE).remove();
}

