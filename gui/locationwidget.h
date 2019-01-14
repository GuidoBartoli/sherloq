#ifndef LOCATIONWIDGET_H
#define LOCATIONWIDGET_H

#include <QtWidgets>
#include "toolwidget.h"
#include "ExifTool.h"

class LocationWidget : public ToolWidget
{
    Q_OBJECT

public:
    LocationWidget(const QString &fileName, QWidget* parent = 0);

private:
    void showError(int type, ExifTool* exif, TagInfo* info);

};

#endif // LOCATIONWIDGET_H
