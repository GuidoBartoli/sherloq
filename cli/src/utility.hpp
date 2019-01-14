#ifndef UTILITY_H
#define UTILITY_H

#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <fstream>
#include <opencv2/opencv.hpp>
#include <boost/filesystem.hpp>

using namespace std;
using namespace cv;

// ----------------------------------------------------------------------------------- //

void printMat(string name, const Mat mat, bool math = false, bool size = false)
{
    if (!math)
    {
        printf("\n");
        if (!name.empty())
            printf("%s (size = %dx%d, type = %d)\n", name.c_str(), mat.rows, mat.cols, mat.type());
        if (size)
            return;
        for (int i = 0; i < mat.rows; ++i)
        {
            printf("    ");
            for (int j = 0; j < mat.cols; ++j)
            {
                if (mat.type() == CV_8U)
                    printf("%4d", mat.at<uchar>(i, j));
                else if (mat.type() == CV_8S)
                    printf("%4d", mat.at<char>(i, j));
                else if (mat.type() == CV_16U)
                    printf("%6d", mat.at<ushort>(i, j));
                else if (mat.type() == CV_16S)
                    printf("%6d", mat.at<short>(i, j));
                else if (mat.type() == CV_32S)
                    printf("%d\t", mat.at<int>(i, j));
                else if (mat.type() == CV_32F)
                    printf("%8.03f", mat.at<float>(i, j));
                else if (mat.type() == CV_64F)
                    printf("%f\t", mat.at<double>(i, j));
            }
            printf("\n");
        }
    }
    else
    {
        if (mat.rows == 1 || mat.cols == 1)
        {
            printf("\nListPlot[{\n");
            Mat temp = mat.reshape(0, 1);
            for (int i = 0; i < temp.cols; ++i)
                printf("%d,%d\n", i, temp.at<uchar>(0, i));
            printf("}]\n");
        }
        else
        {
            printf("\nListContourPlot[{");
            for (int i = 0; i < mat.rows; ++i)
            {
                for (int j = 0; j < mat.cols; ++j)
                {
                    printf("{%d,%d,%f}", i + 1, j + 1, mat.at<float>(i, j));
                    if (i != mat.rows - 1 || j != mat.cols - 1)
                        printf(",");
                }
            }
            printf("}]\n");
        }
    }
}

void saveMat(string name, const Mat mat)
{
    Mat img;
    normalize(mat, img, 0, 255, NORM_MINMAX);
    img.convertTo(img, CV_8U);
    imwrite(name, img);
}

void showMat(string name, const Mat mat, float fx = 1, float fy = 1)
{
    Mat img;
    resize(mat, img, Size(), fx, fy, INTER_NEAREST);
    normalize(img, img, 0, 255, NORM_MINMAX);
    img.convertTo(img, CV_8U);
    imshow(name, img);
}

// ----------------------------------------------------------------------------------- //

// trim from start (in place)
static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
static inline std::string ltrimmed(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
static inline std::string rtrimmed(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
static inline std::string trimmed(std::string s) {
    trim(s);
    return s;
}

// ----------------------------------------------------------------------------------- //

bool isFilePresent(const string file_name)
{
    FILE* file = fopen(trimmed(file_name).c_str(), "r");
    if (!file)
        return false;
    fclose(file);
    return true;
}

bool isFilenameValid(const string file_name)
{
    return (file_name.find(" ") == string::npos);
}

int executeCommand(const string command, bool suppress = true)
{
    string exec = trimmed(command);
    if (suppress)
        exec += " > /dev/null 1> /dev/null 2> /dev/null";
    return system(exec.c_str());
}

string grabOutput(const string command)
{
    char buffer[128];
    std::string result = "";
    std::shared_ptr<FILE> pipe(popen(command.c_str(), "r"), pclose);
    if (!pipe)
        return string("");
    while (!feof(pipe.get())) {
        if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
    }
    return result;
}

int deleteFile(const string file_name)
{
    return executeCommand("rm -frd " + file_name);
}

int deleteFolder(const string folder_name)
{
    return deleteFile(folder_name);
}

int createFolder(const string folder_name)
{
    return executeCommand("mkdir " + folder_name);
}

int copyFile(const string source, const string dest)
{
    return executeCommand("cp " + source + " " + dest);
}

int moveFile(const string source, const string dest)
{
    return executeCommand("mv " + source + " " + dest);
}

bool checkImage(const string image_file)
{
    return (!imread(image_file).empty());
}

string getFileFromPath(string file_path)
{
    return boost::filesystem::path(file_path).filename().string();
}

int getFileSize(const string file_name)
{
    std::ifstream in(file_name.c_str(), std::ifstream::ate | std::ifstream::binary);
    return (in.tellg());
}

void getImageSize(const string file_name, int &height, int &width, float &mpx)
{
    Mat image = imread(file_name);
    height = image.rows;
    width = image.cols;
    mpx = height * width / 1000000.0;
}

string getFileDate(const string file_name)
{
    struct stat attrib;
    stat(file_name.c_str(), &attrib);
    char date[40];
    strftime(date, 40, "%A %d-%m-%Y at %H:%M:%S", localtime(&(attrib.st_mtime)));
    return string(date);
}

string bytesToHuman(int bytes, bool si = true)
{
    const int UNIT = si ? 1000 : 1024;
    if (bytes < UNIT)
        return (to_string(bytes) + " B");
    int exp = (int)(log(bytes) / log(UNIT));
    char pre = string(si ? "kMGTPE" : "KMGTPE").at(exp - 1);
    char human[128];
    sprintf(human, "%.1f %cB", bytes / pow(UNIT, exp), pre);
    return string(human);
}

// ----------------------------------------------------------------------------------- //

const double PI = acos(-1);

Mat buildContrastLUT(uchar low = 0, uchar high = 90)
{
    Mat lut(256, 1, CV_8U);
    if (low > high)
        swap(low, high);
    for (int i = 0; i < 256; ++i)
    {
        if (i <= low)
            lut.at<uchar>(i, 0) = 0;
        else if (i >= high)
            lut.at<uchar>(i, 0) = 255;
        else
            lut.at<uchar>(i, 0) = (uchar)((i - low) * (255.0 / (high - low)));
    }
    return lut;
}

Mat padImage(const Mat image, int block_size)
{
    int right = 0, bottom = 0;
    if (image.cols % block_size != 0)
        right = ((image.cols / block_size) * block_size + block_size) - image.cols;
    if (image.rows % block_size != 0)
        bottom = ((image.rows / block_size) * block_size + block_size) - image.rows;
    if (right == 0 && bottom == 0)
        return image.clone();
    Mat padded;
    if (image.channels() == 1)
        copyMakeBorder(image, padded, 0, bottom, 0, right, BORDER_CONSTANT, Scalar(0));
    else
        copyMakeBorder(image, padded, 0, bottom, 0, right, BORDER_CONSTANT, Scalar(0,0,0));
    return padded;
}

int runExifTool(const string options, const string image_file, const string output_file)
{
    int code;
    string command = "exiftool -m " + options + " " + image_file;
    if (output_file.empty())
        code = executeCommand(command, true);
    else
        code = executeCommand(command + " > " + output_file, false);
    return code;
}

void createPlot(const Mat values, int size, Mat &plot)
{
    plot = Mat(values.cols / 2, values.cols, CV_8U, Scalar(0));
    float v0 = 1;
    for (int i = 0; i < values.cols; ++i)
    {
        float v1 = values.at<float>(0, i);
        Point p0(i - 1, (int)(v0 * plot.rows));
        Point p1(i, (int)(v1 * plot.rows));
        line(plot, p0, p1, Scalar(255));
        v0 = v1;
    }
    resize(plot, plot, Size(), size, size, INTER_NEAREST);
}

#endif // UTILITY_H