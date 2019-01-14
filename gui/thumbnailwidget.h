#ifndef THUMBNAILWIDGET_H
#define THUMBNAILWIDGET_H

#include <QtWidgets>
#include "toolwidget.h"
#include "ExifTool.h"

class ThumbnailWidget : public ToolWidget
{
    Q_OBJECT

public:
    ThumbnailWidget(const QString &fileName, QWidget* parent = Q_NULLPTR);

private:
    void showError(int type, ExifTool* exif, TagInfo* info);
};

#endif // THUMBNAILWIDGET_H
