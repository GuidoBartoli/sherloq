#include "metadatawidget.h"
#include "ExifTool.h"

MetadataWidget::MetadataWidget(const QString &fileName, QWidget* parent) : ToolWidget(parent)
{
    tableWidget = new QTableWidget(0, 2);
    QStringList headers;
    headers << tr("Property") << tr("Value");
    tableWidget->setHorizontalHeaderLabels(headers);
    //tableWidget->setSelectionMode(QAbstractItemView::ExtendedSelection);
    tableWidget->setSelectionMode(QAbstractItemView::SingleSelection);

    treeWidget = new QTreeWidget();
    QStringList labels;
    labels << tr("Group") << tr("Property") << tr("Value");
    treeWidget->setColumnCount(labels.length());
    treeWidget->setHeaderLabels(labels);

    ExifTool* exif = new ExifTool();
    TagInfo* info = exif->ImageInfo(fileName.toStdString().c_str(), "-m");
//    char *out = exif->GetOutput();
//    if (out) qDebug() << out;
//    char *err = exif->GetError();
//    if (err) qDebug() << err;
    if (info)
    {
        int count = 0;
        QStringList groupLevels;
        groupLevels.append(QString());
        groupLevels.append(QString());
        groupLevels.append(QString());
        for (TagInfo* i = info; i; i = i->next)
        {
            QString name = QString::fromUtf8(i->name);
            QString group1 = QString::fromUtf8(i->group[0]);
            QString group2 = QString::fromUtf8(i->group[1]);
            QString group3 = QString::fromUtf8(i->group[2]);
            QString description = QString::fromUtf8(i->desc);
            QString value = QString::fromUtf8(i->value);
            if (name == "SourceFile" || name == "ExifToolVersion" || value.isEmpty() ||
                    value.contains(("binary data"), Qt::CaseInsensitive))
                continue;
            tableWidget->setRowCount(++count);
            tableWidget->setItem(count - 1, 0, new QTableWidgetItem(description));
            tableWidget->setItem(count - 1, 1, new QTableWidgetItem(value));
            QString tooltip = group1 + " --> " + group2 + " --> " + group3 + " --> " + name;
            tableWidget->item(count - 1, 0)->setToolTip(tooltip);
            QFont font = tableWidget->item(count - 1, 0)->font();
            font.setBold(true);
            tableWidget->item(count - 1, 0)->setFont(font);

            QTreeWidgetItem* item1;
            if (group1 != groupLevels.at(0))
            {
                item1 = new QTreeWidgetItem();
                item1->setText(0, group1);
                groupLevels[0] = group1;
                groupLevels[1] = QString();
                groupLevels[2] = QString();
            }
            else
                item1 = treeWidget->topLevelItem(treeWidget->topLevelItemCount() - 1);
            item1->setFont(0, font);
            QTreeWidgetItem* item2;
            if (group2 != groupLevels.at(1))
            {
                item2 = new QTreeWidgetItem(item1);
                item2->setText(0, group2);
                groupLevels[1] = group2;
                groupLevels[2] = QString();
            }
            else
                item2 = item1->child(item1->childCount() - 1);
            item2->setFont(0, font);
            QTreeWidgetItem* item3;
            if (group3 != groupLevels.at(2))
            {
                item3 = new QTreeWidgetItem(item2);
                item3->setText(0, group3);
                groupLevels[2] = group3;
            }
            else
                item3 = item2->child(item2->childCount() - 1);
            item3->setFont(0, font);
            QTreeWidgetItem* item4 = new QTreeWidgetItem(item3);
            item4->setText(1, description);
            item4->setText(2, value);
            treeWidget->addTopLevelItem(item1);
        }
        tableWidget->resizeColumnsToContents();
        tableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);
        tableWidget->setAlternatingRowColors(true);
        treeWidget->expandAll();
        for (int j = 0; j < treeWidget->columnCount(); ++j)
            treeWidget->resizeColumnToContents(j);
        treeWidget->setAnimated(true);
        treeWidget->setAlternatingRowColors(true);
    }
    delete info;
    delete exif;
    connect(tableWidget, &QTableWidget::itemDoubleClicked,
            this,        &MetadataWidget::copyCell);

    tableRadio = new QRadioButton(tr("Table"));
    connect(tableRadio, &QRadioButton::toggled,
            this,       &MetadataWidget::changeView);
    treeRadio = new QRadioButton(tr("Tree"));
    connect(treeRadio, &QRadioButton::toggled,
            this,      &MetadataWidget::changeView);
    expandButton = new QPushButton();
    expandButton->setText(tr("Expand all"));
    connect(expandButton, &QPushButton::clicked,
            treeWidget,   &QTreeWidget::expandAll);
    collapseButton = new QPushButton();
    collapseButton->setText(tr("Collapse all"));
    connect(collapseButton, &QPushButton::clicked,
            treeWidget,     &QTreeWidget::collapseAll);
    QPushButton* exportButton = new QPushButton(tr("Export..."));
    connect(exportButton, &QPushButton::clicked,
            this,         &MetadataWidget::exportData);
    QHBoxLayout* topLayout = new QHBoxLayout();
    topLayout->addWidget(new QLabel(tr("View:")));
    topLayout->addWidget(tableRadio);
    topLayout->addWidget(treeRadio);
    topLayout->addWidget(expandButton);
    topLayout->addWidget(collapseButton);
    topLayout->addWidget(exportButton);
    topLayout->addStretch();

    searchEdit = new QLineEdit();
    connect(searchEdit, &QLineEdit::textChanged,
            this,       &MetadataWidget::checkSearch);
    QToolButton* clearButton = new QToolButton();
    clearButton->setText(tr("Clear"));
    connect(clearButton, &QToolButton::clicked,
            searchEdit,  &QLineEdit::clear);
    QToolButton* prevButton = new QToolButton();
    prevButton->setText(tr("Previous"));
    prevButton->setShortcut(QKeySequence::FindPrevious);
    connect(prevButton, &QToolButton::clicked,
            this,       &MetadataWidget::searchBackward);
    QToolButton* nextButton = new QToolButton();
    nextButton->setText(tr("Next"));
    nextButton->setShortcut(QKeySequence::FindNext);
    connect(nextButton, &QToolButton::clicked,
            this,       &MetadataWidget::searchForward);
    highlightButton = new QToolButton();
    highlightButton->setText(tr("Highlight"));
    highlightButton->setCheckable(true);
    connect(highlightButton, &QToolButton::toggled,
            this,            &MetadataWidget::highlightSearch);
    infoLabel = new QLabel();
    QFont font = infoLabel->font();
    font.setItalic(true);
    infoLabel->setFont(font);

    QHBoxLayout* searchLayout = new QHBoxLayout();
    searchLayout->addWidget(new QLabel(tr("Search:")));
    searchLayout->addWidget(searchEdit);
    searchLayout->addWidget(clearButton);
    searchLayout->addWidget(prevButton);
    searchLayout->addWidget(nextButton);
    searchLayout->addWidget(highlightButton);
    searchLayout->addWidget(infoLabel);
    searchLayout->addStretch();
    searchWidget = new QWidget();
    searchWidget->setLayout(searchLayout);

    QVBoxLayout* mainLayout = new QVBoxLayout();
    mainLayout->addLayout(topLayout);
    mainLayout->addWidget(tableWidget);
    mainLayout->addWidget(treeWidget);
    mainLayout->addWidget(searchWidget);
    setLayout(mainLayout);

    setMinimumSize(600, 300);

    tableRadio->setChecked(true);
}

void MetadataWidget::changeView()
{
    tableWidget->setVisible(tableRadio->isChecked());
    searchWidget->setVisible(tableRadio->isChecked());
    treeWidget->setVisible(treeRadio->isChecked());
    expandButton->setEnabled(treeRadio->isChecked());
    collapseButton->setEnabled(treeRadio->isChecked());
}

void MetadataWidget::checkSearch()
{
    searchFrom(0, 0, +1);
}

void MetadataWidget::exportData()
{
    QString fileName = QFileDialog::getSaveFileName(this, QApplication::applicationName(),
                                                    QString(), tr("Text files (*.txt)"));
    if (fileName.isEmpty())
        return;
    // FIXME: Usare QFileDialog non statico ed impostare setDefaultSuffix("txt")!
    if (QFileInfo(fileName).suffix() != "txt")
        fileName.append(".txt");
    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QMessageBox::critical(this, QApplication::applicationName(),
                              tr("Error while writing to file!"));
        return;
    }

    int maxLength = std::numeric_limits<int>::min();
    for (int i = 0; i < tableWidget->rowCount(); ++i)
    {
        int length = tableWidget->item(i, 0)->text().length();
        if (length > maxLength)
            maxLength = length;
    }

    QTextStream textStream(&file);
    for (int i = 0; i < tableWidget->rowCount(); ++i)
    {
        QString tag = tableWidget->item(i, 0)->text();
        QString value = tableWidget->item(i, 1)->text();
        int numDots = maxLength - tag.length() - 1;
        textStream << tag;
        if (numDots > 0)
            textStream << " ";
        for (int j = 0; j < numDots; ++j)
            textStream << ".";
        textStream << " : " << value << "\r\n";
    }
    file.close();
}

void MetadataWidget::searchFrom(int i0, int j0, int dir)
{
    highlightSearch(highlightButton->isChecked());
    QString pattern = searchEdit->text();
    if (pattern.length() < 2)
    {
        tableWidget->setCurrentCell(-1, -1);
        infoLabel->clear();
        return;
    }

    int i = i0;
    int j = j0;
    while (true)
    {
        while (true)
        {
            QTableWidgetItem* item = tableWidget->item(i, j);
            if (item->text().contains(pattern, Qt::CaseInsensitive))
            {
                tableWidget->setCurrentCell(i, j);
                searchEdit->setStyleSheet("color: #000000");
                return;
            }
            j = dir > 0 ? j + 1 : j - 1;
            if ((dir > 0) && (j == tableWidget->columnCount()))
            {
                j = 0;
                break;
            }
            else if ((dir < 0) && (j == -1))
            {
                j = tableWidget->columnCount() - 1;
                break;
            }
        }
        i = dir > 0 ? i + 1 : i - 1;
        if (((dir > 0) && (i == tableWidget->rowCount())) || ((dir < 0) && (i == -1)))
        {
            searchEdit->setStyleSheet("color: #FF0000");
            infoLabel->setText(tr("%1 of table reached").arg(dir > 0 ? tr("End") : tr("Begin")));
            break;
        }
    }

}

void MetadataWidget::searchForward()
{
    int i0 = tableWidget->currentRow();
    int j0 = tableWidget->currentColumn();
    if (++j0 == tableWidget->columnCount())
    {
        ++i0;
        j0 = 0;
    }
    if (i0 == tableWidget->rowCount())
        return;
    searchFrom(i0, j0, +1);
}

void MetadataWidget::searchBackward()
{
    int i0 = tableWidget->currentRow();
    int j0 = tableWidget->currentColumn();
    if (--j0 == -1)
    {
        --i0;
        j0 = tableWidget->columnCount() - 1;
    }
    if (i0 == -1)
        return;
    searchFrom(i0, j0, -1);
}

void MetadataWidget::highlightSearch(bool enabled)
{
    QString pattern = searchEdit->text();
    if (pattern.isEmpty())
        return;
    int count = 0;
    for (int i = 0; i < tableWidget->rowCount(); ++i)
        for (int j = 0; j < tableWidget->columnCount(); ++j)
            if (enabled && pattern.length() >= 2 &&
                    tableWidget->item(i, j)->text().contains(pattern, Qt::CaseInsensitive))
            {
                tableWidget->item(i, j)->setBackground(Qt::yellow);
                ++count;
            }
            else
                tableWidget->item(i, j)->setBackground(Qt::transparent);
    if (enabled)
        infoLabel->setText(tr("%1 occurence(s) found").arg(count).arg(pattern));
    else
        infoLabel->clear();
}

void MetadataWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->matches(QKeySequence::Find))
        searchEdit->setFocus();
    else
        ToolWidget::keyPressEvent(event);
}

void MetadataWidget::copyCell(QTableWidgetItem *item)
{
    QApplication::clipboard()->setText(item->text());
    emit messageToShow(tr("[FILE::Metadata] Cell contents copied to clipboard"));
}
