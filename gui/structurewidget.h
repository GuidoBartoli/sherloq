#ifndef STRUCTUREWIDGET_H
#define STRUCTUREWIDGET_H

#include <QtWidgets>
#include "toolwidget.h"

class StructureWidget : public ToolWidget
{

    Q_OBJECT

    const QString TEMP_FILE = QDir::tempPath() + "/structure.html";

public:
    StructureWidget(const QString &fileName, QWidget* parent = Q_NULLPTR);
    ~StructureWidget();

};

#endif // STRUCTUREWIDGET_H
