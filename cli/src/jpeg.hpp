#ifndef JPEG_HPP
#define JPEG_HPP

#include <opencv2/opencv.hpp>
#include <dirent.h>
#include "utility.hpp"

using namespace std;
using namespace cv;

const int DCT_SIZE = 8;
const int TABLE_SIZE = DCT_SIZE * DCT_SIZE;
const int DCT_OFFSET = 128;

const Mat ZIG_ZAG = (Mat_<int>(TABLE_SIZE, 2) <<
        0, 0,   0, 1,   1, 0,   2, 0,   1, 1,   0, 2,   0, 3,   1, 2,
        2, 1,   3, 0,   4, 0,   3, 1,   2, 2,   1, 3,   0, 4,   0, 5,
        1, 4,   2, 3,   3, 2,   4, 1,   5, 0,   6, 0,   5, 1,   4, 2,
        3, 3,   2, 4,   1, 5,   0, 6,   0, 7,   1, 6,   2, 5,   3, 4,
        4, 3,   5, 2,   6, 1,   7, 0,   7, 1,   6, 2,   5, 3,   4, 4,
        3, 5,   2, 6,   1, 7,   2, 7,   3, 6,   4, 5,   5, 4,   6, 3,
        7, 2,   7, 3,   6, 4,   5, 5,   4, 6,   3, 7,   4, 7,   5, 6,
        6, 5,   7, 4,   7, 5,   6, 6,   5, 7,   6, 7,   7, 6,   7, 7);

void getQTs(int quality, Mat &luma_qt, Mat &chroma_qt, bool baseline = true)
{
    const Mat LUMA_QT = (Mat_<double>(DCT_SIZE, DCT_SIZE) <<
        16,  11,  10,  16,  24,  40,  51,  61,
        12,  12,  14,  19,  26,  58,  60,  55,
        14,  13,  16,  24,  40,  57,  69,  56,
        14,  17,  22,  29,  51,  87,  80,  62,
        18,  22,  37,  56,  68, 109, 103,  77,
        24,  35,  55,  64,  81, 104, 113,  92,
        49,  64,  78,  87, 103, 121, 120, 101,
        72,  92,  95,  98, 112, 100, 103,  99);
    const Mat CHROMA_QT = (Mat_<double>(DCT_SIZE, DCT_SIZE) <<
        17,  18,  24,  47,  99,  99,  99,  99,
        18,  21,  26,  66,  99,  99,  99,  99,
        24,  26,  56,  99,  99,  99,  99,  99,
        47,  66,  99,  99,  99,  99,  99,  99,
        99,  99,  99,  99,  99,  99,  99,  99,
        99,  99,  99,  99,  99,  99,  99,  99,
        99,  99,  99,  99,  99,  99,  99,  99,
        99,  99,  99,  99,  99,  99,  99,  99);

    if (quality <= 0)
        quality = 1;
    else if (quality > 100)
        quality = 100;
    if (quality < 50)
        quality = 5000 / quality;
    else
        quality = 200 - quality * 2;

    luma_qt = (LUMA_QT * quality + 50) / 100;
    luma_qt.convertTo(luma_qt, baseline ? CV_8U : CV_16U);
    luma_qt.setTo(1, luma_qt <= 0);
    chroma_qt = (CHROMA_QT * quality + 50) / 100;
    chroma_qt.convertTo(chroma_qt, baseline ? CV_8U : CV_16U);
    chroma_qt.setTo(1, chroma_qt <= 0);
}

bool findNextMarker(FILE* file, uchar byte1, uchar byte2, uchar byte3)
{
    while (!feof(file))
        if (fgetc(file) == byte1 && fgetc(file) == byte2 && fgetc(file) == byte3)
            return true;
    return false;
}

int qualityEstimation(const string filename, Mat &luma, float &luma_q, Mat &chroma, float &chroma_q, int &jpeg_quality)
{
    const string TEMP_FILE = "temp.jpg";

    const uchar JPEG_MRK = 0xFF;
    const uchar JPEG_SOI = 0xD8;
    const uchar JPEG_DQT = 0xDB;
    const uchar JPEG_PAD = 0x00;

    const int MAX_TABLES = 2;
    const int LEN_OFFSET = 2;
    const int START_QUAL = 0;
    const int END_QUAL   = 100;
    const int QUAL_STEP  = 1;
    const int NUM_STEPS  = (END_QUAL - START_QUAL) / QUAL_STEP + 1;
    const int LUMA_IDX   = 0;
    const int CHROMA_IDX = 1;

    FILE* file = fopen(filename.c_str(), "rb");
    int first = fgetc(file);
    int second = fgetc(file);
    fclose(file);
    if (first != JPEG_MRK || second != JPEG_SOI)
        return -1;
    copyFile(filename, TEMP_FILE);
    if (runExifTool("-all=", TEMP_FILE, "") != 0)
        return -2;
    deleteFile(TEMP_FILE + "_original");

    file = fopen(TEMP_FILE.c_str(), "rb");
    bool found = false;
    luma.create(DCT_SIZE, DCT_SIZE, CV_8U);
    chroma.create(DCT_SIZE, DCT_SIZE, CV_8U);
    while (!feof(file))
    {
        if (!findNextMarker(file, JPEG_MRK, JPEG_DQT, JPEG_PAD))
            break;
        int length = fgetc(file) - LEN_OFFSET;
        if ((length % (TABLE_SIZE + 1)) != 0 || length <= 0)
            continue;
        while (length > 0)
        {
            int type = fgetc(file);
            --length;
            int index = type & 0x0f;
            if (index >= MAX_TABLES)
                break;
            for (int k = 0; k < TABLE_SIZE; ++k)
            {
                uchar c = fgetc(file);
                --length;
                if (feof(file))
                    break;
                int i = ZIG_ZAG.at<int>(k, 0);
                int j = ZIG_ZAG.at<int>(k, 1);
                if (index == LUMA_IDX)
                    luma.at<uchar>(i, j) = c;
                else if (index == CHROMA_IDX)
                    chroma.at<uchar>(i, j) = c;
            }
            if (luma.at<uchar>(0, 0) != 0)
                found = true;
        }
    }
    fclose(file);
    deleteFile(TEMP_FILE);
    if (!found)
        return -3;

    float luma_avg = mean(luma)[0];
    luma_avg *= TABLE_SIZE;
    luma_avg -= luma.at<uchar>(0, 0);
    luma_avg /= TABLE_SIZE - 1;
    luma_q = (1 - luma_avg / 255) * 100;
    float chroma_avg = mean(chroma)[0];
    chroma_avg *= TABLE_SIZE;
    chroma_avg -= chroma.at<uchar>(0, 0);
    chroma_avg /= TABLE_SIZE - 1;
    chroma_q = (1 - chroma_avg / 255) * 100;
    Mat search(1, NUM_STEPS, CV_32F, Scalar(0));
    int count = 0;
    for (int q = START_QUAL; q <= END_QUAL; q += QUAL_STEP)
    {
        Mat luma0;
        Mat chroma0;
        getQTs(q, luma0, chroma0);
        float diff = mean(abs(luma - luma0))[0] + mean(abs(chroma - chroma0))[0];
        search.at<float>(0, count++) = diff;
    }
    double min_val;
    Point min_loc;
    minMaxLoc(search, &min_val, 0, &min_loc);
    jpeg_quality = START_QUAL + QUAL_STEP * min_loc.x;

    return 0;
}

// ----------------------------------------------------------------------------------- //

Mat jpegCompress(const Mat input, int quality, bool color = true)
{
    vector<int> params;
    params.push_back(CV_IMWRITE_JPEG_QUALITY);
    params.push_back(quality);
    vector<uchar> buffer;
    imencode(".jpg", input, buffer, params);
    return (imdecode(buffer, color ? CV_LOAD_IMAGE_COLOR : CV_LOAD_IMAGE_GRAYSCALE));
}

void errorLevelAnalysis(const Mat input, Mat &output)
{
    const int QUALITY = 75;
    const float SCALE = 20;

    Mat compressed = jpegCompress(input, QUALITY);
    compressed.convertTo(compressed, CV_32F);
    Mat original;
    input.convertTo(original, CV_32F);
    output = abs(compressed - original) * SCALE;
    output.convertTo(output, CV_8U);
}

// ----------------------------------------------------------------------------------- //

int getFirstDigit(int number)
{
    char buffer[12];
    sprintf(buffer, "%d", abs(number));
    return (buffer[0] - '0');
}

void firstDigitFeatures(const Mat input, Mat &features)
{
    const int SKIP_DC    = 1;
    const int NUM_MODES  = 20;
    const int NUM_DIGITS = 10;
    const vector<int> DIGITS = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    const Mat SCALING = (Mat_<float>(NUM_MODES * DIGITS.size(), 2) <<
        0.023303, 0.513002,    0.005404, 0.306991,    0.004053, 0.465721,
        0.003377, 0.346842,    0.000675, 0.119216,    0.000675, 0.364404,
        0.000338, 0.416751,    0.000000, 0.309017,    0.000000, 0.363728,
        0.031746, 0.538332,    0.008781, 0.313745,    0.004728, 0.170550,
        0.003715, 0.321851,    0.000000, 0.298210,    0.002364, 0.340763,
        0.001689, 0.448159,    0.000000, 0.348531,    0.000338, 0.183384,
        0.028369, 0.557582,    0.001351, 0.295846,    0.000000, 0.315096,
        0.000252, 0.340892,    0.000338, 0.133739,    0.000000, 0.337048,
        0.000338, 0.348450,    0.000000, 0.339075,    0.000000, 0.236744,
        0.009456, 0.520770,    0.000338, 0.331982,    0.000675, 0.166160,
        0.000338, 0.340426,    0.000000, 0.299223,    0.000000, 0.350895,
        0.000338, 0.437352,    0.000000, 0.368794,    0.000000, 0.076663,
        0.027693, 0.499831,    0.003040, 0.340763,    0.000000, 0.374536,
        0.000338, 0.338399,    0.000338, 0.336035,    0.000338, 0.491726,
        0.000000, 0.394461,    0.000000, 0.353597,    0.000675, 0.455589,
        0.026342, 0.537318,    0.000000, 0.291118,    0.000000, 0.363728,
        0.000000, 0.170550,    0.000000, 0.400878,    0.000000, 0.360014,
        0.000000, 0.088146,    0.000000, 0.349544,    0.000000, 0.284363,
        0.028369, 0.554205,    0.000675, 0.299223,    0.000000, 0.355623,
        0.000000, 0.344478,    0.000000, 0.106045,    0.000000, 0.343465,
        0.000000, 0.346505,    0.000000, 0.399865,    0.000000, 0.282675,
        0.018237, 0.532590,    0.004390, 0.302938,    0.000000, 0.340763,
        0.000000, 0.342452,    0.000000, 0.389058,    0.000000, 0.271868,
        0.000000, 0.416751,    0.000000, 0.460655,    0.000000, 0.370483,
        0.040865, 0.605539,    0.000675, 0.295846,    0.000000, 0.344433,
        0.000000, 0.357987,    0.000000, 0.122256,    0.000000, 0.361027,
        0.000000, 0.373522,    0.000000, 0.379264,    0.000000, 0.251942,
        0.039176, 0.612631,    0.000000, 0.348193,    0.000000, 0.297872,
        0.000000, 0.357649,    0.000000, 0.345154,    0.000000, 0.102668,
        0.000000, 0.388720,    0.000000, 0.313745,    0.000000, 0.404255,
        0.021952, 0.561297,    0.000000, 0.334347,    0.000000, 0.355285,
        0.000000, 0.153664,    0.000000, 0.460993,    0.000000, 0.273556,
        0.000000, 0.380952,    0.000000, 0.088821,    0.000000, 0.387707,
        0.020601, 0.559608,    0.000000, 0.292131,    0.000000, 0.356636,
        0.000000, 0.181695,    0.000000, 0.475177,    0.000000, 0.359338,
        0.000000, 0.093549,    0.000000, 0.373185,    0.000000, 0.291456,
        0.016548, 0.535630,    0.000000, 0.322188,    0.000338, 0.289092,
        0.000000, 0.366093,    0.000000, 0.147923,    0.000000, 0.479905,
        0.000000, 0.109085,    0.000000, 0.377913,    0.000000, 0.073286,
        0.006417, 0.544748,    0.000000, 0.541709,    0.000000, 0.231341,
        0.000000, 0.246538,    0.000000, 0.350718,    0.000000, 0.137116,
        0.000000, 0.449510,    0.000000, 0.095576,    0.000000, 0.083418,
        0.019926, 0.590341,    0.000000, 0.467072,    0.000000, 0.441405,
        0.000000, 0.405268,    0.000000, 0.162445,    0.000000, 0.127660,
        0.000000, 0.216481,    0.000000, 0.376224,    0.000000, 0.068558,
        0.009119, 0.563323,    0.000000, 0.558595,    0.000000, 0.315772,
        0.000000, 0.237082,    0.000000, 0.366591,    0.000000, 0.121243,
        0.000000, 0.090510,    0.000000, 0.476528,    0.000000, 0.069909,
        0.008443, 0.576494,    0.000000, 0.562648,    0.000000, 0.231003,
        0.000000, 0.262749,    0.000000, 0.365755,    0.000000, 0.110773,
        0.000000, 0.466059,    0.000000, 0.061128,    0.000000, 0.104357,
        0.008781, 0.547112,    0.000000, 0.545424,    0.000000, 0.208713,
        0.000000, 0.377575,    0.000000, 0.165485,    0.000000, 0.111787,
        0.000000, 0.474164,    0.000000, 0.318811,    0.000000, 0.365755,
        0.006417, 0.548126,    0.000000, 0.565012,    0.000000, 0.193178,
        0.000000, 0.373522,    0.000000, 0.103006,    0.000000, 0.088821,
        0.000000, 0.471125,    0.000000, 0.308342,    0.000000, 0.373522,
        0.021952, 0.555893,    0.000000, 0.582236,    0.000000, 0.230665,
        0.000000, 0.210402,    0.000000, 0.369363,    0.000000, 0.119554,
        0.000000, 0.453901,    0.000000, 0.052009,    0.000000, 0.113813);

    Mat gray = padImage(input, DCT_SIZE);
    cvtColor(gray, gray, CV_BGR2GRAY);
    gray.convertTo(gray, CV_32F);
    Mat digits_pdf(NUM_MODES, NUM_DIGITS, CV_32F);
    digits_pdf.setTo(0);
    for (int i = 0; i < input.rows - DCT_SIZE; i += DCT_SIZE)
    {
        for (int j = 0; j < input.cols - DCT_SIZE; j += DCT_SIZE)
        {
            Mat block(gray, Rect(j, i, DCT_SIZE, DCT_SIZE));
            block -= 128;
            Mat coeffs(DCT_SIZE, DCT_SIZE, CV_32F);
            dct(block, coeffs);
            coeffs.convertTo(coeffs, CV_16S);
            for (int m = 0; m < NUM_MODES; ++m)
            {
                int i0 = ZIG_ZAG.at<int>(m + SKIP_DC, 0);
                int j0 = ZIG_ZAG.at<int>(m + SKIP_DC, 1);
                short coeff = coeffs.at<short>(i0, j0);
                int d = getFirstDigit(coeff);
                ++digits_pdf.at<float>(m, d);
            }
        }
    }

    features.create(1, NUM_MODES * DIGITS.size(), CV_32F);
    int f = 0;
    for (int m = 0; m < NUM_MODES; ++m)
    {
        float digits_sum = 0;
        for (int d = 0; d < NUM_DIGITS; ++d)
            digits_sum += digits_pdf.at<float>(m, d);
        for (int d = 0; d < NUM_DIGITS; ++d)
        {
            if (find(DIGITS.begin(), DIGITS.end(), d) != DIGITS.end())
            {
                float pdf = digits_pdf.at<float>(m, d) / digits_sum;
                float min_val = SCALING.at<float>(f, 0);
                float max_val = SCALING.at<float>(f, 1);
                features.at<float>(0, f) = 2 * (pdf - min_val) / (max_val - min_val) - 1;
                ++f;
            }
        }
    }
}

int doubleCompression(const Mat image, float &probability)
{
    const string MODEL_FILE = string("../svm/first_digit.yml");
    // const float MIN_DIST = -5.73780;
    const float MIN_DIST = -1.76819;
    const float MAX_DIST = +1.36294;
    // const float MAX_DIST = +3.19925;

    if (!isFilePresent(MODEL_FILE))
        return -1;

    CvSVM svm;
    svm.load(MODEL_FILE.c_str());
    Mat features;
    firstDigitFeatures(image, features);
    float dist = svm.predict(features, true);
    if (dist > 0)
    {
        if (dist > MAX_DIST)
            probability = -100;
        else
            probability = -dist / MAX_DIST * 100;
    }
    else
    {
        if (dist < MIN_DIST)
            probability = +100;
        else
            probability = dist / MIN_DIST * 100;
    }
    return 0;
}

// ----------------------------------------------------------------------------------- //

void compressionGhosts(const Mat input, int &n, int &l1, float &p1, int &l2, float &p2, int &q1, int &q2, Mat &plot)
{
    const int MIN_QUALITY  = 3;
    const int MAX_QUALITY  = 97;
    const int QUALITY_STEP = 1;
    const int KERNEL_SIZE  = 5;
    const int PLOT_RESIZE  = 10;

    const float MIN_PEAK = 0.218;
    const float THR_PEAK = 0.521;
    const float MAX_PEAK = 0.823;

    const bool DEBUG = false;

    Mat gray;
    cvtColor(input, gray, CV_BGR2GRAY);
    Mat original;
    gray.convertTo(original, CV_32F);

    vector<int> jpeg_qualities;
    for (int q = MAX_QUALITY; q >= MIN_QUALITY; q -= QUALITY_STEP)
        jpeg_qualities.push_back(q);
    Mat jpeg_losses(1, jpeg_qualities.size(), CV_32F);
    for (uint i = 0; i < jpeg_qualities.size(); ++i)
    {
        Mat compressed = jpegCompress(input, jpeg_qualities.at(i), false);
        compressed.convertTo(compressed, CV_32F);
        Scalar loss = mean(abs(original - compressed));
        jpeg_losses.at<float>(0, i) = loss[0];
    }
    normalize(jpeg_losses, jpeg_losses, 0, 1, NORM_MINMAX);
    createPlot(jpeg_losses, PLOT_RESIZE, plot);

    Mat peak_response;
    Mat peak_kernel = (Mat_<float>(1, KERNEL_SIZE) << 0, 0, -1, 0, 0);
    matchTemplate(jpeg_losses, peak_kernel, peak_response, CV_TM_CCOEFF_NORMED);
    Mat padding(1, KERNEL_SIZE / 2, CV_32F, Scalar(0));
    hconcat(padding, peak_response, peak_response);
    hconcat(peak_response, padding, peak_response);
    double peak_val;
    Point peak_loc;
    minMaxLoc(peak_response, 0, &peak_val, 0, &peak_loc);

    // Scalar mean, stddev;
    // meanStdDev(peak_response, mean, stddev);
    // float dev = 3*stddev[0];
    // p1 = max_peak;
    // return;

    if (DEBUG)
    {
        printf("\nListLinePlot[{{");
        for (uint i = 0; i < jpeg_qualities.size(); ++i)
            printf("{%d,%.04f}, ", jpeg_qualities.at(i), jpeg_losses.at<float>(0, i));
        printf("},\n{");
        for (uint i = 0; i < jpeg_qualities.size(); ++i)
            printf("{%d,%.04f}, ", jpeg_qualities.at(i), peak_response.at<float>(0, i));
        printf("}},PlotRange->All]\n");
    }

    if (peak_val < THR_PEAK)
    {
        l1 = l2 = p2 = -1;
        p1 = (peak_val - THR_PEAK) / (MIN_PEAK - THR_PEAK) * 100;
        n = 0;
    }
    else
    {
        l1 = jpeg_qualities.at(peak_loc.x);
        p1 = (peak_val - THR_PEAK) / (MAX_PEAK - THR_PEAK) * 100;
        peak_response.at<float>(0, peak_loc.x) = peak_response.at<float>(0, peak_loc.x - 1) = peak_response.at<float>(0, peak_loc.x + 1) = -1;
        minMaxLoc(peak_response, 0, &peak_val, 0, &peak_loc);
        if (peak_val > THR_PEAK)
        {
            l2 = jpeg_qualities.at(peak_loc.x);
            p2 = (peak_val - THR_PEAK) / (MAX_PEAK - THR_PEAK) * 100;
            if (p1 > p2)
            {
                swap(l1, l2);
                swap(p1, p2);
            }
            n = 2;
            for (int i = 0; i < peak_response.cols; ++i)
                if (peak_response.at<float>(0, i) > THR_PEAK)
                ++n;
        }
        else
        {
            l2 = p2 = -1;
            n = 1;
        }
    }
    if (p1 > 100)
        p1 = 100;
    if (p2 > 100)
        p2 = 100;
    q1 = MIN_QUALITY + KERNEL_SIZE / 2;
    q2 = MAX_QUALITY - KERNEL_SIZE / 2;
}

// ----------------------------------------------------------------------------------- //

#endif // JPEG_HPP