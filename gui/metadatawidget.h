#ifndef METADATAWIDGET_H
#define METADATAWIDGET_H

#include <QtWidgets>
#include "toolwidget.h"

class MetadataWidget : public ToolWidget
{
    Q_OBJECT

public:
    MetadataWidget(const QString &fileName, QWidget* parent = Q_NULLPTR);

protected:
    void keyPressEvent(QKeyEvent* event) Q_DECL_OVERRIDE;

private:
    QRadioButton* tableRadio;
    QRadioButton* treeRadio;
    QPushButton*  expandButton;
    QPushButton*  collapseButton;
    QTableWidget* tableWidget;
    QTreeWidget*  treeWidget;
    QLineEdit*    searchEdit;
    QWidget*      searchWidget;
    QLabel*       infoLabel;
    QToolButton*  highlightButton;
    void searchFrom(int i0, int j0, int dir);

private slots:
    void changeView();
    void searchForward();
    void searchBackward();
    void highlightSearch(bool enabled);
    void checkSearch();
    void exportData();
    void copyCell(QTableWidgetItem* item);

};

#endif // METADATAWIDGET_H
