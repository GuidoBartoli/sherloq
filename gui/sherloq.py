import sys

import cv2 as cv
from PySide2.QtCore import Qt, QSettings, QFileInfo
from PySide2.QtGui import QKeySequence, QIcon
from PySide2.QtWidgets import (
    QApplication,
    QMainWindow,
    QMdiArea,
    QMdiSubWindow,
    QDockWidget,
    QAction,
    QMessageBox,
    QFileDialog)

from adjust import AdjustWidget
from digest import DigestWidget
from reverse import ReverseWidget
from echo import EchoWidget
from editor import EditorWidget
from ela import ElaWidget
from fourier import FourierWidget
from gradient import GradientWidget
from location import LocationWidget
from metadata import MetadataWidget
from minmax import MinMaxWidget
from multiple import MultipleWidget
from noise import NoiseWidget
from original import OriginalWidget
from pca import PcaWidget
from planes import PlanesWidget
from quality import QualityWidget
from space import SpaceWidget
from stats import StatsWidget
from structure import StructureWidget
from thumbnail import ThumbWidget
from tools import ToolTree
from utility import modify_font


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        QApplication.setApplicationName('Sherloq')
        QApplication.setOrganizationName('Guido Bartoli')
        QApplication.setOrganizationDomain('www.guidobartoli.com')
        QApplication.setApplicationVersion(ToolTree().version)
        QApplication.setWindowIcon(QIcon('icons/sherloq_white.png'))
        self.setWindowTitle('{} {}'.format(QApplication.applicationName(), QApplication.applicationVersion()))
        self.mdi_area = QMdiArea()
        self.setCentralWidget(self.mdi_area)
        self.filename = None
        self.image = None
        modify_font(self.statusBar(), bold=True)

        tree_dock = QDockWidget(self.tr('TOOLS'), self)
        tree_dock.setObjectName('tree_dock')
        tree_dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)
        self.addDockWidget(Qt.LeftDockWidgetArea, tree_dock)
        self.tree_widget = ToolTree()
        self.tree_widget.setObjectName('tree_widget')
        self.tree_widget.itemDoubleClicked.connect(self.open_tool)
        tree_dock.setWidget(self.tree_widget)

        tree_action = tree_dock.toggleViewAction()
        tree_action.setToolTip(self.tr('Toggle toolset visibility'))
        tree_action.setText(self.tr('Show tools'))
        tree_action.setShortcut(QKeySequence(Qt.Key_Tab))
        tree_action.setObjectName('tree_action')
        tree_action.setIcon(QIcon('icons/tools.svg'))

        load_action = QAction(self.tr('&Load image...'), self)
        load_action.setToolTip(self.tr('Choose an image to analyze'))
        load_action.setShortcut(QKeySequence.Open)
        load_action.triggered.connect(self.load_file)
        load_action.setObjectName('load_action')
        load_action.setIcon(QIcon('icons/load.svg'))

        quit_action = QAction(self.tr('&Quit'), self)
        quit_action.setToolTip(self.tr('Exit from the program'))
        quit_action.setShortcut(QKeySequence.Quit)
        quit_action.triggered.connect(self.close)
        quit_action.setObjectName('quit_action')
        quit_action.setIcon(QIcon('icons/quit.svg'))

        tabbed_action = QAction(self.tr('&Tabbed'), self)
        tabbed_action.setToolTip(self.tr('Toggle tabbed view for window area'))
        tabbed_action.setShortcut(QKeySequence(Qt.Key_F10))
        tabbed_action.setCheckable(True)
        tabbed_action.triggered.connect(self.toggle_view)
        tabbed_action.setObjectName('tabbed_action')
        tabbed_action.setIcon(QIcon('icons/tabbed.svg'))

        prev_action = QAction(self.tr('&Previous'), self)
        prev_action.setToolTip(self.tr('Select the previous tool window'))
        prev_action.setShortcut(QKeySequence.PreviousChild)
        prev_action.triggered.connect(self.mdi_area.activatePreviousSubWindow)
        prev_action.setObjectName('prev_action')
        prev_action.setIcon(QIcon('icons/previous.svg'))

        next_action = QAction(self.tr('&Next'), self)
        next_action.setToolTip(self.tr('Select the next tool window'))
        next_action.setShortcut(QKeySequence.NextChild)
        next_action.triggered.connect(self.mdi_area.activateNextSubWindow)
        next_action.setObjectName('next_action')
        next_action.setIcon(QIcon('icons/next.svg'))

        tile_action = QAction(self.tr('&Tile'), self)
        tile_action.setToolTip(self.tr('Arrange windows into non-overlapping views'))
        tile_action.setShortcut(QKeySequence(Qt.Key_F11))
        tile_action.triggered.connect(self.mdi_area.tileSubWindows)
        tile_action.setObjectName('tile_action')
        tile_action.setIcon(QIcon('icons/tile.svg'))

        cascade_action = QAction(self.tr('&Cascade'), self)
        cascade_action.setToolTip(self.tr('Arrange windows into overlapping views'))
        cascade_action.setShortcut(QKeySequence(Qt.Key_F12))
        cascade_action.triggered.connect(self.mdi_area.cascadeSubWindows)
        cascade_action.setObjectName('cascade_action')
        cascade_action.setIcon(QIcon('icons/cascade.svg'))

        close_action = QAction(self.tr('Close &All'), self)
        close_action.setToolTip(self.tr('Close all open tool windows'))
        close_action.setShortcut(QKeySequence(Qt.CTRL + Qt.SHIFT + Qt.Key_W))
        close_action.triggered.connect(self.mdi_area.closeAllSubWindows)
        close_action.setObjectName('close_action')
        close_action.setIcon(QIcon('icons/close.svg'))

        full_action = QAction(self.tr('Full screen'), self)
        full_action.setToolTip(self.tr('Enter full screen mode'))
        full_action.setShortcut(QKeySequence.FullScreen)
        full_action.triggered.connect(self.showFullScreen)
        full_action.setObjectName('full_action')
        full_action.setIcon(QIcon('icons/full.svg'))

        normal_action = QAction(self.tr('Normal view'), self)
        full_action.setToolTip(self.tr('Switch to normal screen mode'))
        normal_action.setShortcut(QKeySequence(Qt.CTRL + Qt.Key_F12))
        normal_action.triggered.connect(self.showNormal)
        normal_action.setObjectName('normal_action')
        normal_action.setIcon(QIcon('icons/normal.svg'))

        about_action = QAction(self.tr('&About...'), self)
        about_action.setToolTip(self.tr('Display informations about this program'))
        about_action.setShortcut(QKeySequence.HelpContents)
        about_action.triggered.connect(self.show_about)
        about_action.setObjectName('about_action')
        about_action.setIcon(QIcon('icons/sherloq_alpha.png'))

        about_qt_action = QAction(self.tr('About &Qt'), self)
        about_qt_action.setToolTip(self.tr('Display informations about the Qt Framework'))
        about_qt_action.triggered.connect(QApplication.aboutQt)
        about_qt_action.setIcon(QIcon('icons/Qt.svg'))

        file_menu = self.menuBar().addMenu(self.tr('File'))
        file_menu.addAction(load_action)
        file_menu.addAction(quit_action)

        window_menu = self.menuBar().addMenu(self.tr('&Window'))
        window_menu.addAction(tree_action)
        window_menu.addAction(tabbed_action)
        window_menu.addSeparator()
        window_menu.addAction(tile_action)
        window_menu.addAction(cascade_action)
        window_menu.addSeparator()
        window_menu.addAction(prev_action)
        window_menu.addAction(next_action)
        window_menu.addAction(close_action)
        window_menu.addSeparator()
        window_menu.addAction(full_action)
        window_menu.addAction(normal_action)

        help_menu = self.menuBar().addMenu(self.tr('&Help'))
        help_menu.addAction(about_action)
        help_menu.addAction(about_qt_action)

        main_toolbar = self.addToolBar(self.tr('&Toolbar'))
        main_toolbar.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        main_toolbar.addAction(load_action)
        main_toolbar.addAction(tree_action)
        main_toolbar.addSeparator()
        main_toolbar.addAction(prev_action)
        main_toolbar.addAction(next_action)
        main_toolbar.addSeparator()
        main_toolbar.addAction(tile_action)
        main_toolbar.addAction(cascade_action)
        main_toolbar.addAction(tabbed_action)
        main_toolbar.addAction(close_action)
        # main_toolbar.addSeparator()
        # main_toolbar.addAction(normal_action)
        # main_toolbar.addAction(full_action)
        main_toolbar.setObjectName('main_toolbar')

        settings = QSettings()
        settings.beginGroup('main_window')
        self.restoreGeometry(settings.value('geometry'))
        self.restoreState(settings.value('state'))
        settings.endGroup()

        prev_action.setEnabled(False)
        next_action.setEnabled(False)
        tile_action.setEnabled(False)
        cascade_action.setEnabled(False)
        close_action.setEnabled(False)
        tabbed_action.setEnabled(False)
        self.tree_widget.setEnabled(False)
        self.show_message(self.tr('Ready'))

    def closeEvent(self, event):
        settings = QSettings()
        settings.beginGroup('main_window')
        settings.setValue('geometry', self.saveGeometry())
        settings.setValue('state', self.saveState())
        settings.endGroup()
        super(MainWindow, self).closeEvent(event)

    def load_file(self):
        settings = QSettings()
        self.filename = QFileDialog.getOpenFileName(
            self, self.tr('Load image'), settings.value('load_folder'),
            self.tr('Supported formats (*.jpg *.jpeg *.jpe *.jp2 *.png *.tif *.tiff, *.bmp, *.gif)'))[0]
        if not self.filename:
            return
        image = cv.imread(self.filename, cv.IMREAD_COLOR)
        if image is None:
            QMessageBox.critical(self, self.tr('Error'), self.tr('Unable to load image!'))
            return
        if image.shape[2] > 3:
            QMessageBox.warning(self, self.tr('Warning'), self.tr('Embedded alpha channel discarded'))
            image = cv.cvtColor(image, cv.COLOR_BGRA2BGR)
        self.image = image

        path = QFileInfo(self.filename).absolutePath()
        name = QFileInfo(self.filename).fileName()
        settings.setValue('load_folder', path)
        self.findChild(ToolTree, 'tree_widget').setEnabled(True)
        self.findChild(QAction, 'prev_action').setEnabled(True)
        self.findChild(QAction, 'next_action').setEnabled(True)
        self.findChild(QAction, 'tile_action').setEnabled(True)
        self.findChild(QAction, 'cascade_action').setEnabled(True)
        self.findChild(QAction, 'close_action').setEnabled(True)
        self.findChild(QAction, 'tabbed_action').setEnabled(True)
        self.setWindowTitle('[{}] - {} {}'.format(
            name, QApplication.applicationName(), QApplication.applicationVersion()))
        self.show_message(self.tr('Image "{}" successfully loaded'.format(name)))

        # FIXME: disable_bold della chiusura viene chiamato DOPO open_tool e nell'albero la voce NON diventa neretto
        self.mdi_area.closeAllSubWindows()
        self.open_tool(self.tree_widget.topLevelItem(0).child(0), None)

    def open_tool(self, item, _):
        if not item.data(0, Qt.UserRole):
            return
        group = item.data(0, Qt.UserRole + 1)
        tool = item.data(0, Qt.UserRole + 2)
        for sub_window in self.mdi_area.subWindowList():
            if sub_window.windowTitle() == item.text(0):
                sub_window.setWindowState(Qt.WindowActive)
                sub_window.setFocus()
                return

        tool_widget = None
        if group == 0:
            if tool == 0:
                tool_widget = OriginalWidget(self.image)
            elif tool == 1:
                tool_widget = DigestWidget(self.filename, self.image)
            elif tool == 2:
                tool_widget = EditorWidget()
            elif tool == 3:
                tool_widget = ReverseWidget()
            else:
                return
        elif group == 1:
            if tool == 0:
                tool_widget = StructureWidget(self.filename)
            elif tool == 1:
                tool_widget = MetadataWidget(self.filename)
            elif tool == 2:
                tool_widget = ThumbWidget(self.filename, self.image)
            elif tool == 3:
                tool_widget = LocationWidget(self.filename)
            else:
                return
        elif group == 2:
            if tool == 2:
                tool_widget = AdjustWidget(self.image)
            elif tool == 3:
                tool_widget = FourierWidget(self.image)
            else:
                return
        elif group == 3:
            if tool == 0:
                tool_widget = QualityWidget(self.filename)
            elif tool == 1:
                tool_widget = ElaWidget(self.image)
            elif tool == 2:
                tool_widget = MultipleWidget(self.image)
            else:
                return
        elif group == 4:
            if tool == 1:
                tool_widget = PcaWidget(self.image)
            elif tool == 2:
                tool_widget = StatsWidget(self.image)
            elif tool == 3:
                tool_widget = SpaceWidget(self.image)
            else:
                return
        elif group == 5:
            if tool == 0:
                tool_widget = GradientWidget(self.image)
            elif tool == 1:
                tool_widget = EchoWidget(self.image)
            else:
                return
        elif group == 6:
            if tool == 0:
                tool_widget = NoiseWidget(self.image)
            elif tool == 1:
                tool_widget = MinMaxWidget(self.image)
            elif tool == 2:
                tool_widget = PlanesWidget(self.image)
        else:
            return
        # FIXME: Aggiungere un metodo init e dopo fare il connect, sennò i messaggi del costruttore non si vedono
        tool_widget.info_message.connect(self.show_message)
        tool_widget.help_clicked.connect(self.show_help)

        sub_window = QMdiSubWindow()
        sub_window.setWidget(tool_widget)
        sub_window.setWindowTitle(item.text(0))
        sub_window.setObjectName(item.text(0))
        sub_window.setAttribute(Qt.WA_DeleteOnClose)
        sub_window.setWindowIcon(QIcon('icons/{}.svg'.format(group)))
        self.mdi_area.addSubWindow(sub_window)
        sub_window.show()
        sub_window.destroyed.connect(self.disable_bold)
        self.tree_widget.set_bold(item.text(0), enabled=True)

    def disable_bold(self, item):
        self.tree_widget.set_bold(item.windowTitle(), enabled=False)

    def toggle_view(self, tabbed):
        if tabbed:
            self.mdi_area.setViewMode(QMdiArea.TabbedView)
            self.mdi_area.setTabsClosable(True)
            self.mdi_area.setTabsMovable(True)
        else:
            self.mdi_area.setViewMode(QMdiArea.SubWindowView)
        self.findChild(QAction, 'tile_action').setEnabled(not tabbed)
        self.findChild(QAction, 'cascade_action').setEnabled(not tabbed)

    def show_about(self):
        message = '<h2>{} {}</h2>'.format(QApplication.applicationName(), QApplication.applicationVersion())
        message += '<h3>A digital image forensic toolkit</h3>'
        message += '<p>author: <a href="{}">{}</a></p>'.format(
            QApplication.organizationDomain(), QApplication.organizationName())
        message += '<p>source: <a href="https://github.com/GuidoBartoli/sherloq">GitHub repository</a></p>'
        QMessageBox.about(self, self.tr('About'), message)

    def show_message(self, message):
        self.statusBar().showMessage(message, 10000)

    def show_help(self, tool):
        self.show_message('help from {}'.format(tool))
        pass


if __name__ == '__main__':
    application = QApplication(sys.argv)
    mainwindow = MainWindow()
    mainwindow.show()
    sys.exit(application.exec_())
