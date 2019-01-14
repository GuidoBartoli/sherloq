#include "mainwindow.h"
#include "cvmatandqimage.h"

#include "originalwidget.h"
#include "digestwidget.h"
#include "metadatawidget.h"
#include "structurewidget.h"
#include "thumbnailwidget.h"
#include "locationwidget.h"
#include "qualitywidget.h"
#include "ghostwidget.h"
#include "gradientwidget.h"
#include "sweepwidget.h"
#include "rgbpixelwidget.h"
#include "minmaxwidget.h"
#include "spacewidget.h"
#include "elawidget.h"
#include "colorpcawidget.h"
#include "adjustmentwidget.h"
#include "separationwidget.h"
#include "magnifierwidget.h"
#include "doublewidget.h"
#include "comparisonwidget.h"
#include "echowidget.h"
#include "contrastwidget.h"
#include "snrwidget.h"
#include "histogramwidget.h"

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
{
    QApplication::setApplicationName("Sherlock");
    QApplication::setOrganizationName("EsseGi_Soft");
    QApplication::setOrganizationDomain("www.essegisoft.com");
    QApplication::setApplicationVersion("0.22a");
    setWindowTitle(QApplication::applicationName());

    mdiArea = new QMdiArea();
    setCentralWidget(mdiArea);
    //QFont font = statusBar()->font();
    //font.setBold(true);
    //statusBar()->setFont(font);
    move(QApplication::desktop()->screen()->rect().center() - rect().center());

    QDockWidget* treeDock = new QDockWidget(tr("TOOLS"), this);
    treeDock->setObjectName("ToolDock");
    treeDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    toolTree = new ToolTree();
    treeDock->setWidget(toolTree);
    addDockWidget(Qt::LeftDockWidgetArea, treeDock);
    connect(toolTree, &ToolTree::itemDoubleClicked,
            this,     &MainWindow::openTool);
    QAction* treeAction = treeDock->toggleViewAction();
    treeAction->setToolTip(tr("Toggle toolset visibility"));
    treeAction->setText(tr("Tools"));
    treeAction->setShortcut(Qt::Key_Tab);

    QAction* loadAction = new QAction(tr("&Load image..."), this);
    loadAction->setToolTip(tr("Choose an image to analyze"));
    loadAction->setShortcut(QKeySequence::Open);
    connect(loadAction, &QAction::triggered,
            this,       &MainWindow::loadFile);

    QAction* quitAction = new QAction(tr("&Quit"), this);
    quitAction->setToolTip(tr("Exit from the program"));
    quitAction->setShortcut(QKeySequence::Quit);
    connect(quitAction, &QAction::triggered,
            this,       &MainWindow::close);

    QAction* tabbedAction = new QAction(tr("&Tabbed view"), this);
    tabbedAction->setToolTip(tr("Toggle tabbed view for tool windows"));
    tabbedAction->setShortcut(Qt::Key_F10);
    tabbedAction->setCheckable(true);
    connect(tabbedAction, &QAction::toggled,
            this,         &MainWindow::toggleView);

    tileAction = new QAction(tr("&Tile"), this);
    tileAction->setToolTip(tr("Arrange windows into non-overlapping views"));
    tileAction->setShortcut(Qt::Key_F11);
    connect(tileAction, &QAction::triggered,
            mdiArea,    &QMdiArea::tileSubWindows);

    cascadeAction = new QAction(tr("&Cascade"), this);
    cascadeAction->setToolTip(tr("Arrange windows into overlapping views"));
    cascadeAction->setShortcut(Qt::Key_F12);
    connect(cascadeAction, &QAction::triggered,
            mdiArea,       &QMdiArea::cascadeSubWindows);

    nextAction = new QAction(tr("&Next"), this);
    nextAction->setToolTip(tr("Select the next tool window"));
    nextAction->setShortcut(QKeySequence::NextChild);
    connect(nextAction, &QAction::triggered,
            mdiArea,    &QMdiArea::activateNextSubWindow);

    previousAction = new QAction(tr("&Previous"), this);
    previousAction->setToolTip(tr("Select the previous tool window"));
    previousAction->setShortcut(QKeySequence::PreviousChild);
    connect(previousAction, &QAction::triggered,
            mdiArea,        &QMdiArea::activatePreviousSubWindow);

    closeAllAction = new QAction(tr("Close &All"), this);
    closeAllAction->setToolTip(tr("Close all open tool windows"));
    closeAllAction->setShortcut(Qt::CTRL + Qt::SHIFT + Qt::Key_W);
    connect(closeAllAction, &QAction::triggered,
            mdiArea,        &QMdiArea::closeAllSubWindows);

    QAction* aboutAction = new QAction(tr("&About..."), this);
    aboutAction->setToolTip(tr("Display informations about this program"));
    aboutAction->setShortcut(QKeySequence::HelpContents);
    connect(aboutAction, &QAction::triggered,
            this,        &MainWindow::showAbout);

    QAction* aboutQtAction = new QAction(tr("About &Qt"), this);
    aboutQtAction->setToolTip(tr("Display informations about the Qt Framework"));
    connect(aboutQtAction, &QAction::triggered,
            qApp,          &QApplication::aboutQt);

    QMenu* fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(loadAction);
    fileMenu->addAction(quitAction);
    QMenu* windowMenu = menuBar()->addMenu(tr("&Window"));
    windowMenu->addAction(treeAction);
    windowMenu->addAction(tabbedAction);
    windowMenu->addSeparator();
    windowMenu->addAction(tileAction);
    windowMenu->addAction(cascadeAction);
    windowMenu->addSeparator();
    windowMenu->addAction(nextAction);
    windowMenu->addAction(previousAction);
    windowMenu->addAction(closeAllAction);
    QMenu* helpMenu = menuBar()->addMenu(tr("&Help"));
    helpMenu->addAction(aboutAction);
    helpMenu->addAction(aboutQtAction);

    QToolBar* mainToolBar = addToolBar(tr("&Toolbar"));
    mainToolBar->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    mainToolBar->addAction(loadAction);
    mainToolBar->addSeparator();
    mainToolBar->addAction(treeAction);
    mainToolBar->addSeparator();
    mainToolBar->addAction(tileAction);
    mainToolBar->addAction(cascadeAction);
    mainToolBar->addSeparator();
    mainToolBar->addAction(nextAction);
    mainToolBar->addAction(previousAction);
    mainToolBar->setObjectName("Toolbar");

    QSettings settings;
    settings.beginGroup("mainwindow");
    restoreGeometry(settings.value("geometry").toByteArray());
    restoreState(settings.value("state").toByteArray());
    settings.endGroup();

    toolTree->setEnabled(false);
    previousAction->setEnabled(false);
    nextAction->setEnabled(false);
    tileAction->setEnabled(false);
    cascadeAction->setEnabled(false);
    closeAllAction->setEnabled(false);

    showMessage(tr("%1 engine started").arg(QApplication::applicationName()));
}

MainWindow::~MainWindow()
{
    QSettings settings;
    settings.beginGroup("mainwindow");
    settings.setValue("geometry", saveGeometry());
    settings.setValue("state", saveState());
    settings.endGroup();
}

void MainWindow::openTool(QTreeWidgetItem* item, int)
{
    if (!item->data(0, Qt::UserRole).toBool())
        return;
    QString toolName = item->text(0);
    for (int i = 0; i < mdiArea->subWindowList().length(); ++i)
    {
        if (mdiArea->subWindowList().at(i)->windowTitle() == toolName)
        {
            mdiArea->subWindowList().at(i)->setFocus();
            return;
        }
    }
    QFont font = item->font(0);
    font.setBold(true);
    item->setFont(0, font);

    int cat = item->data(0, Qt::UserRole + 1).toInt();
    int tool = item->data(0, Qt::UserRole + 2).toInt();
    ToolWidget* toolWidget;
    if (cat == 0) // General
    {
        if (tool == 0) // Original Image
            toolWidget = new OriginalWidget(fileName, imageCv, this);
        else if (tool == 1) // Image Digest
            toolWidget = new DigestWidget(fileName, this);
        else
            return;
    }
    else if (cat == 1)
    {
        if (tool == 0)
            toolWidget = new MetadataWidget(fileName, this);
        else if (tool == 1)
            toolWidget = new StructureWidget(fileName, this);
        else if (tool == 2)
            toolWidget = new ThumbnailWidget(fileName, this);
        else if (tool == 3)
            toolWidget = new LocationWidget(fileName, this);
        else
            return;
    }
    else if (cat == 2)
    {
        if (tool == 0)
            toolWidget = new MagnifierWidget(imageCv, this);
        else if (tool == 1)
            toolWidget = new AdjustmentWidget(imageCv, this);
        else if (tool == 2)
            toolWidget = new SweepWidget(imageCv, this);
        else if (tool == 3)
            toolWidget = new ComparisonWidget(imageCv, this);
        else
            return;
    }
    else if (cat == 3)
    {
        if (tool == 0)
            toolWidget = new QualityWidget(fileName, this);
        else if (tool == 1)
            toolWidget = new GhostWidget(imageCv, this);
        else if (tool == 2)
            toolWidget = new DoubleWidget(imageCv, this);
        else if (tool == 3)
            toolWidget = new ElaWidget(imageCv, this);
        else
            return;
    }
    else if (cat == 4)
    {
        if (tool == 0)
            toolWidget = new HistogramWidget(imageCv, this);
        else if (tool == 1)
            toolWidget = new SpaceWidget(imageCv, this);
        else if (tool == 2)
            toolWidget = new ColorPcaWidget(imageCv, this);
        else if (tool == 3)
            toolWidget = new RgbPixelWidget(imageCv, this);
        else
            return;
    }
    else if (cat == 5)
    {
        if (tool == 0)
            toolWidget = new GradientWidget(imageCv, this);
        else if (tool == 2)
            toolWidget = new EchoWidget(imageCv, this);
        else
            return;
    }
    else if (cat == 6)
    {
        if (tool == 0)
            toolWidget = new SeparationWidget(imageCv, this);
        else if (tool == 1)
            toolWidget = new MinMaxWidget(imageCv, this);
        else if (tool == 2)
            toolWidget = new SnrWidget(imageCv, this);
        else
            return;
    }
    else if (cat == 7)
    {
        if (tool == 0)
            toolWidget = new ContrastWidget(imageCv, this);
        else
            return;
    }
    else
        return;

    QMdiSubWindow* subWindow = new QMdiSubWindow();
    subWindow->setWidget(toolWidget);
    subWindow->setWindowTitle(toolName);
    subWindow->setAttribute(Qt::WA_DeleteOnClose);
    mdiArea->addSubWindow(subWindow);
    subWindow->show();

    connect(toolWidget, &ToolWidget::messageToShow,
            this,       &MainWindow::showMessage);
    connect(subWindow, &QMdiSubWindow::destroyed,
            this,      &MainWindow::disableBold);
}

void MainWindow::loadFile()
{
    QSettings settings;
    QString tempFile = QFileDialog::getOpenFileName(
                this, tr("Load image"),
                settings.value("last_folder").toString(),
                tr("Supported images (*.jpg *.jpeg *.png *.tif *.tiff)"));
    if (tempFile.isEmpty())
        return;
    QString path = QFileInfo(tempFile).path();
    settings.setValue("last_folder", path);

    cv::Mat tempImage = cv::imread(tempFile.toStdString(), CV_LOAD_IMAGE_COLOR);
    if (tempImage.empty())
    {
        QMessageBox::critical(this, QApplication::applicationName(), tr("Unable to load image!"));
        return;
    }
    //qDebug() << tempImage.type() << tempImage.channels();
    if (tempImage.channels() > 3)
    {
        QMessageBox::warning(this, QApplication::applicationName(), tr("Embedded alpha channel discarded"));
        tempImage.convertTo(tempImage, CV_BGRA2BGR);
    }

    mdiArea->closeAllSubWindows();
    imageCv = tempImage;
    fileName = tempFile;

    toolTree->setEnabled(true);
    previousAction->setEnabled(true);
    nextAction->setEnabled(true);
    tileAction->setEnabled(true);
    cascadeAction->setEnabled(true);
    closeAllAction->setEnabled(true);

    showMessage(tr("Image \"%1\" successfully loaded").arg(QFileInfo(fileName).fileName()));

    openTool(toolTree->getItem(0, 0), 0);
}

void MainWindow::showAbout()
{
    QMessageBox::
            about(this, tr("About"), tr("<h2>%1 %2").
                  arg(QApplication::applicationName()).
                  arg(QApplication::applicationVersion()) + "</h2>" +
                  tr("<h3>A digital image forensic tool</h3>"
                     "<p>Copyright 2017, <a href='%1'>%2</a>. "
                     "All rights reserved.</p>").
                  arg(QApplication::organizationDomain()).
                  arg(QApplication::organizationName()));
}

void MainWindow::showMessage(const QString &message)
{
    statusBar()->showMessage(message, 5000);
}

void MainWindow::disableBold(QObject* object)
{
    QString toolName = static_cast<QMdiSubWindow*>(object)->windowTitle();
    toolTree->disableBold(toolName);
}

void MainWindow::toggleView(bool tabbed)
{
    if (tabbed)
    {
        mdiArea->setViewMode(QMdiArea::TabbedView);
        mdiArea->setTabsClosable(true);
        mdiArea->setTabsMovable(true);
    }
    else
        mdiArea->setViewMode(QMdiArea::SubWindowView);
}
