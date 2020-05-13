from PySide2.QtCore import Qt, QSettings, QFileInfo
from PySide2.QtGui import QKeySequence
from PySide2.QtWidgets import (
    QApplication,
    QAbstractItemView,
    QVBoxLayout,
    QFileDialog,
    QTableWidgetItem,
    QTableWidget,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QToolButton)
from tabulate import tabulate

from pyexiftool import exiftool
from utility import modify_font
from widget import ToolWidget


class MetadataWidget(ToolWidget):
    def __init__(self, filename, parent=None):
        super(MetadataWidget, self).__init__(parent)
        self.table_widget = QTableWidget(0, 3)
        headers = [self.tr('Group'), self.tr('Property'), self.tr('Value')]
        self.table_widget.setHorizontalHeaderLabels(headers)
        self.table_widget.setSelectionMode(QAbstractItemView.SingleSelection)
        rows = 0
        last = None
        with exiftool.ExifTool('pyexiftool/exiftool/exiftool') as et:
            metadata = et.get_metadata(filename)
            for tag, value in metadata.items():
                ignore = ['SourceFile', 'ExifTool:ExifTool', 'File:FileName', 'File:Directory',
                          'File:FileSize', 'File:FileModifyDate', 'File:FileAccessDate', 'File:FileInodeChangeDate',
                          'File:FileType', 'File:FilePermissions', 'File:FileTypeExtension', 'File:MIMEType']
                if not value or any(t in tag for t in ignore):
                    continue
                value = str(value).replace(', use -b option to extract', '')
                group, desc = tag.split(':')
                self.table_widget.setRowCount(rows + 1)
                if last is None or group != last:
                    self.table_widget.setItem(rows, 0, QTableWidgetItem(group))
                    last = group
                self.table_widget.setItem(rows, 1, QTableWidgetItem(desc))
                self.table_widget.setItem(rows, 2, QTableWidgetItem(value))
                modify_font(self.table_widget.item(rows, 0), bold=True)
                rows += 1
        self.table_widget.resizeColumnsToContents()
        self.table_widget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table_widget.setAlternatingRowColors(True)
        self.table_widget.itemDoubleClicked.connect(self.copy)

        self.search_edit = QLineEdit()
        self.search_edit.textChanged.connect(self.check)
        clear_button = QToolButton()
        clear_button.setText(self.tr('Clear'))
        clear_button.setShortcut(QKeySequence.DeleteCompleteLine)
        clear_button.clicked.connect(self.search_edit.clear)
        prev_button = QToolButton()
        prev_button.setText(self.tr('Previous'))
        prev_button.setShortcut(QKeySequence.FindPrevious)
        prev_button.clicked.connect(self.backward)
        next_button = QToolButton()
        next_button.setText(self.tr('Next'))
        prev_button.setShortcut(QKeySequence.FindNext)
        next_button.clicked.connect(self.forward)
        high_button = QToolButton()
        high_button.setText(self.tr('Highlight'))
        high_button.setCheckable(True)
        high_button.toggled.connect(self.highlight)
        self.info_label = QLabel()
        modify_font(self.info_label, italic=True)
        export_button = QToolButton()
        export_button.setText(self.tr('Export...'))
        modify_font(export_button, bold=True)
        export_button.clicked.connect(self.export)

        search_layout = QHBoxLayout()
        search_layout.addWidget(QLabel(self.tr('Search:')))
        search_layout.addWidget(self.search_edit)
        search_layout.addWidget(self.info_label)
        search_layout.addWidget(clear_button)
        search_layout.addWidget(prev_button)
        search_layout.addWidget(next_button)
        search_layout.addWidget(high_button)
        search_layout.addStretch()
        search_layout.addWidget(export_button)

        layout = QVBoxLayout()
        layout.addWidget(self.table_widget)
        layout.addLayout(search_layout)
        self.setLayout(layout)
        self.setMinimumSize(800, 600)

    def check(self):
        self.search(0, 0, 1)

    def forward(self):
        i0 = self.table_widget.currentRow()
        j0 = self.table_widget.currentColumn() + 1
        if j0 == self.table_widget.columnCount():
            i0 += 1
            j0 = 0
        if i0 == self.table_widget.rowCount():
            return
        self.search(i0, j0, +1)

    def backward(self):
        i0 = self.table_widget.currentRow()
        j0 = self.table_widget.currentColumn() - 1
        if j0 == -1:
            i0 -= 1
            j0 = self.table_widget.columnCount() - 1
        if i0 == -1:
            return
        self.search(i0, j0, -1)

    def search(self, i0, j0, direction):
        pattern = self.search_edit.text().lower()
        if len(pattern) < 2:
            self.table_widget.setCurrentCell(-1, -1)
            self.info_label.clear()
            return
        i = i0
        j = j0
        rows = self.table_widget.rowCount()
        cols = self.table_widget.columnCount()
        while True:
            while True:
                item = self.table_widget.item(i, j)
                if item is not None and pattern in item.text().lower():
                    self.table_widget.setCurrentCell(i, j)
                    self.search_edit.setStyleSheet('color: #000000')
                    return
                j = j + 1 if direction > 0 else j - 1
                if direction > 0 and j == cols:
                    j = 0
                    break
                if direction < 0 and j == -1:
                    j = cols - 1
                    break
            i = i + 1 if direction > 0 else i - 1
            if (direction > 0 and i == rows) or (direction < 0 and i == -1):
                self.search_edit.setStyleSheet('color: #FF0000')
                self.info_label.setText('{} of table reached'.format('End' if direction > 0 else 'Top'))
                break

    def highlight(self, enabled):
        pattern = self.search_edit.text().lower()
        if not pattern:
            return
        count = 0
        for i in range(self.table_widget.rowCount()):
            for j in range(self.table_widget.columnCount()):
                item = self.table_widget.item(i, j)
                if item is None:
                    continue
                if enabled and len(pattern) >= 2 and pattern in item.text().lower():
                    self.table_widget.item(i, j).setBackground(Qt.yellow)
                    count += 1
                else:
                    self.table_widget.item(i, j).setBackground(Qt.transparent)
        if enabled:
            self.info_label.setText(self.tr('{} occurence(s) found'.format(count)))
        else:
            self.info_label.clear()

    def export(self):
        settings = QSettings()
        filename = QFileDialog.getSaveFileName(
            self, self.tr('Export metadata'), settings.value('save_folder'), self.tr('Text files (*.txt)'))[0]
        if not filename:
            return
        if not filename.endswith('.txt'):
            filename += '.txt'

        rows = self.table_widget.rowCount()
        cols = self.table_widget.columnCount()
        table = [[None for _ in range(cols)] for __ in range(rows)]
        for i in range(rows):
            for j in range(cols):
                item = self.table_widget.item(i, j)
                if item is not None:
                    table[i][j] = item.text()
        with open(filename, 'w') as file:
            file.write(tabulate(table, headers=['Group', 'Property', 'Value'], tablefmt='psql'))
        settings.setValue('save_folder', QFileInfo(filename).absolutePath())

    def copy(self, item):
        QApplication.clipboard().setText(item.text())
        self.message_to_show.emit(self.tr('Cell contents copied to clipboard'))
