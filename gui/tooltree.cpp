#include "tooltree.h"

ToolTree::ToolTree(QWidget *parent) : QTreeWidget(parent)
{
    QStringList categoryNames;
    QList<QStringList> toolNames;
    QList<QStringList> toolInfos;
    QStringList names;
    QStringList infos;

    categoryNames.append(tr("General"));
    names << tr("Original Image") << tr("Image Digest") <<
             tr("Similarity Search") << tr("Automatic Tagging");
    toolNames.append(names);
    names.clear();
    infos.append(tr("Display the unaltered reference image for visual inspection"));
    infos.append(tr("Compute byte and perceptual hashes together with extension ballistics"));
    infos.append(tr("Use online reverse search services for finding similar images on the web"));
    infos.append(tr("Apply deep learning algorithms for automatic picture tagging"));
    toolInfos.append(infos);
    infos.clear();

    categoryNames.append(tr("File"));
    names << tr("Metadata Extraction") << tr("EXIF Structure") <<
             tr("Thumbnail Analysis") << tr("Geolocation Data");
    toolNames.append(names);
    names.clear();
    infos.append(tr("Scan through file metadata and gather all available information"));
    infos.append(tr("Dump the physical EXIF structure and display and interactive view"));
    infos.append(tr("Extract embedded thumbnail (if present) and compare with original"));
    infos.append(tr("Get geo-location data (if present) and display on a world map view"));
    toolInfos.append(infos);
    infos.clear();

    categoryNames.append(tr("Inspection"));
    names << tr("Enhancing Magnifier") << tr("Image Adjustments") <<
             tr("Tonal Range Sweep") << tr("Reference Comparison");
    toolNames.append(names);
    names.clear();
    infos.append(tr("Apply various visual enhancement for better identifying forgeries"));
    infos.append(tr("Apply typical adjustments (contrast, brightness, hue, saturation, ...)"));
    infos.append(tr("Compress tonality over various levels to detect composite artifacts"));
    infos.append(tr("Open a synchronized double view to compare two different pictures"));
    toolInfos.append(infos);
    infos.clear();

    categoryNames.append(tr("JPEG"));
    names << tr("Quality Estimation") << tr("Compression Ghosts") <<
             tr("Double Compression") << tr("Error Level Analysis");
    toolNames.append(names);
    names.clear();
    infos.append(tr("Extract quantization/huffman tables and estimate last saved JPEG quality"));
    infos.append(tr("Use error residuals to detect multiple compressions at different levels"));
    infos.append(tr("Exploit First Digit Statistics to discover potential double compression"));
    infos.append(tr("Identify areas with different compression levels against a fixed quality"));
    toolInfos.append(infos);
    infos.clear();

    categoryNames.append(tr("Colors"));
    names << tr("RGB/HSV 3D Plots") << tr("Color Space Conversion") <<
             tr("Principal Component Analysis") << tr("RGB Pixel Statistics");
    toolNames.append(names);
    names.clear();
    infos.append(tr("Display interactive 2D/3D plots of RGB and HSV pixel data"));
    infos.append(tr("Convert color channels into RGB/HSV/YCbCr/Lab/CMYK color spaces"));
    infos.append(tr("Use color PCA to project RGB values onto a different vector space"));
    infos.append(tr("Compute Minimum/Maximum/Average RGB values for every pixel"));
    toolInfos.append(infos);
    infos.clear();

    categoryNames.append(tr("Luminance"));
    names << tr("Luminance Gradient") << tr("Luminance Detail") <<
             tr("Echo Edge Filter") << tr("Wavelet Analysis");
    toolNames.append(names);
    names.clear();
    infos.append(tr("Analyze brightness variations along X and Y axis of the image"));
    infos.append(tr("Estimate the high-frequency component of the luminance channel"));
    infos.append(tr("Use laplacian filter to reveal artificial out-of-focus zones"));
    infos.append(tr("Use wavelet reconstruction with various coefficient thresholds"));
    toolInfos.append(infos);
    infos.clear();

    categoryNames.append(tr("Noise"));
    names << tr("Noise Separation") << tr("Min/Max Deviation") <<
             tr("SNR Consistency") << tr("Noise Segmentation");
    toolNames.append(names);
    names.clear();
    infos.append(tr("Estimate and separate natural noise component of the image"));
    infos.append(tr("Highlight pixels deviating from block-based Min/Max statistics"));
    infos.append(tr("Evaluate uniformity of signal-to-noise ratio across the image"));
    infos.append(tr("Cluster noise into uniform zones for easier composite detection"));
    toolInfos.append(infos);
    infos.clear();

    categoryNames.append(tr("Tampering"));
    names << tr("Contrast Enhancement") << tr("Clone Detection") <<
             tr("Resampling Detection") << tr("Splicing Detection");
    toolNames.append(names);
    names.clear();
    infos.append(tr("Analyze image histogram for contrast enhancement detection"));
    infos.append(tr("Use feature descriptors for copy/rotate clone area detection"));
    infos.append(tr("Analyze 2D pixle intepolation for detecting resampling traces"));
    infos.append(tr("Use DCT statistics for automatic splicing zone detection"));
    toolInfos.append(infos);
    infos.clear();

    for (int i = 0; i < categoryNames.length(); ++i)
    {
        QTreeWidgetItem* categoryItem = new QTreeWidgetItem();
        categoryItem->setText(0, categoryNames.at(i));
        QFont font = categoryItem->font(0);
        font.setBold(true);
        categoryItem->setFont(0, font);
        categoryItem->setData(0, Qt::UserRole, false);
        for (int j = 0; j < toolNames.at(i).length(); ++j)
        {
            QTreeWidgetItem* toolItem = new QTreeWidgetItem(categoryItem);
            toolItem->setText(0, toolNames.at(i).at(j));
            toolItem->setData(0, Qt::UserRole, true);
            toolItem->setData(0, Qt::UserRole + 1, i);
            toolItem->setData(0, Qt::UserRole + 2, j);
            toolItem->setToolTip(0, toolInfos.at(i).at(j));
        }
        addTopLevelItem(categoryItem);
    }
    expandAll();
    setColumnCount(1);
    header()->setVisible(false);
    setMaximumWidth(300);
}

void ToolTree::disableBold(QString toolName)
{
    for (int i = 0; i < topLevelItemCount(); ++i)
    {
        QTreeWidgetItem* category = topLevelItem(i);
        for (int j = 0; j < category->childCount(); ++j)
        {
            QTreeWidgetItem* item = category->child(j);
            if (item->text(0) == toolName)
            {
                QFont font = item->font(0);
                font.setBold(false);
                item->setFont(0, font);
                return;
            }
        }
    }
}

void ToolTree::chooseTool(int i, int j)
{
    if (i < 0 || i >= topLevelItemCount() || j < 0 || j >= topLevelItem(i)->childCount())
        return;
    /*
    QTreeWidgetItem* item = topLevelItem(i)->child(j);
    QFont font = item->font(0);
    font.setBold(true);
    item->setFont(0, font);
    */
    qDebug() << topLevelItem(i)->child(j)->text(0);
    emit itemDoubleClicked(topLevelItem(i)->child(j), 0);
}

QTreeWidgetItem *ToolTree::getItem(int i, int j)
{
    return topLevelItem(i)->child(j);
}

/*
 *
QString ToolSet::getCategoryName(int index) const
{
    return categoryNames.at(index);
}

int ToolSet::getCategoryIndex(QString name) const
{
    return categoryNames.indexOf(QRegExp(name));
}

QStringList ToolSet::getCategoryNames() const
{
    return categoryNames;
}

QStringList ToolSet::getToolNames(int index) const
{
    return toolNames.at(index);
}

int ToolSet::getCategoryCount() const
{
    return categoryNames.length();
}

int ToolSet::getToolCount(int index) const
{
    return toolNames.at(index).length();
}

QString ToolSet::getToolName(int index1, int index2) const
{
    return toolNames.at(index1).at(index2);
}
*/
