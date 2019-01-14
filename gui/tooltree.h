#ifndef TOOLTREE_H
#define TOOLTREE_H

#include <QtWidgets>

class ToolTree : public QTreeWidget
{
    Q_OBJECT

public:
    explicit ToolTree(QWidget* parent = Q_NULLPTR);
    void disableBold(QString toolName);
    void chooseTool(int i, int j);
    QTreeWidgetItem* getItem(int i, int j);

    /*
    QStringList getCategoyNames()                   const;
    QStringList getToolNames(int index)             const;
    QString     getCategoryName(int index)          const;
    int         getCategoryIndex(QString name)      const;
    QString     getToolName(int index1, int index2) const;
    int         getCategoryCount()                  const;
    int         getToolCount(int index)             const;
    */
};

#endif // TOOLTREE_H
