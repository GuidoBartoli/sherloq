#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtWidgets>
#include <opencv2/opencv.hpp>
#include "tooltree.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget* parent = Q_NULLPTR);
    ~MainWindow();

private:
    QMdiArea* mdiArea;
    ToolTree* toolTree;
    QString   fileName;
    cv::Mat   imageCv;

    QAction* tileAction;
    QAction* cascadeAction;
    QAction* nextAction;
    QAction* previousAction;
    QAction* closeAllAction;

private slots:
    void openTool(QTreeWidgetItem* item, int);
    void loadFile();
    void showAbout();
    void showMessage(const QString &message);
//    void disableBold(QString name);
    void disableBold(QObject* object);
    void toggleView(bool tabbed);

};

#endif // MAINWINDOW_H
