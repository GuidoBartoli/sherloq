#ifndef TEST_HPP
#define TEST_HPP

#include <opencv2/opencv.hpp>
#include <dirent.h>
#include "utility.hpp"
#include "jpeg.hpp"

void bacmFeatures(const Mat input, Mat &features)
{
    return;

    const int  BLOCK_SIZE   = 8;
    const int  MAX_VALUE    = 255 * 2;
    const int  NUM_FEATURES = 14;
    const bool DEBUG        = false;

    Mat gray;
    cvtColor(input, gray, CV_BGR2GRAY);
    float bacm[BLOCK_SIZE][BLOCK_SIZE];
    for (int m = 0; m < BLOCK_SIZE; ++m)
        for (int n = 0; n < BLOCK_SIZE; ++n)
            bacm[m][n] = 0;

    // int step = ((gray.rows * gray.cols) / (BLOCK_SIZE * BLOCK_SIZE)) / 10;
    // int block_count = 0;
    for (int i0 = 0; i0 < gray.rows - BLOCK_SIZE; i0 += BLOCK_SIZE)
    {
        for (int j0 = 0; j0 < gray.cols - BLOCK_SIZE; j0 += BLOCK_SIZE)
        {
            // if (block_count++ % step == 0)
            // {
            //     printf(".");
            //     fflush(stdout);
            // }
            int h1[BLOCK_SIZE][BLOCK_SIZE][MAX_VALUE];
            int h2[BLOCK_SIZE][BLOCK_SIZE][MAX_VALUE];
            for (int m = 0; m < BLOCK_SIZE; ++m)
                for (int n = 0; n < BLOCK_SIZE; ++n)
                    for (int v = 0; v < MAX_VALUE; ++v)
                        h1[m][n][v] = h2[m][n][v] = 0;

            for (int i = i0; i < i0 + BLOCK_SIZE; ++i)
            {
                for (int j = j0; j < j0 + BLOCK_SIZE; ++j)
                {
                    int a = gray.at<uchar>(i, j);
                    int b = gray.at<uchar>(i, j + 1);
                    int c = gray.at<uchar>(i + 1, j);
                    int d = gray.at<uchar>(i + 1, j + 1);
                    int e = gray.at<uchar>(i + 4, j + 4);
                    int f = gray.at<uchar>(i + 4, j + 5);
                    int g = gray.at<uchar>(i + 5, j + 4);
                    int h = gray.at<uchar>(i + 5, j + 5);
                    h1[i - i0][j - j0][abs(a - b - c + d)]++;
                    h2[i - i0][j - j0][abs(e - f - g + h)]++;
                }
            }

            for (int m = 0; m < BLOCK_SIZE; ++m)
                for (int n = 0; n < BLOCK_SIZE; ++n)
                    for (int v = 0; v < MAX_VALUE; ++v)
                        bacm[m][n] += fabs(h1[m][n][v] - h2[m][n][v]) / MAX_VALUE;
        }
    }

    float min_val = numeric_limits<float>::max();
    float max_val = numeric_limits<float>::min();
    for (int m = 0; m < BLOCK_SIZE; ++m)
    {
        for (int n = 0; n < BLOCK_SIZE; ++n)
        {
            if (bacm[m][n] > max_val)
                max_val = bacm[m][n];
            if (bacm[m][n] < min_val)
                min_val = bacm[m][n];
        }
    }
    for (int m = 0; m < BLOCK_SIZE; ++m)
        for (int n = 0; n < BLOCK_SIZE; ++n)
            bacm[m][n] = (bacm[m][n] - min_val) / (max_val - min_val);

    //   01234567
    // 0    |   *
    // 1 R1 V R2*
    // 2    |   *
    // 3 -H-C-H-*
    // 4    |   *
    // 5 R4 V R3*
    // 6    |   *
    // 7 ********

    // H-V central simmetry
    float h1h2c = (fabs(bacm[3][0] - bacm[3][6]) + fabs(bacm[3][1] - bacm[3][5]) + fabs(bacm[3][2] - bacm[3][4]));
    h1h2c = (h1h2c / 3); // - 0.5) * 2;
    float v1v2c = (fabs(bacm[0][3] - bacm[6][3]) + fabs(bacm[1][3] - bacm[5][3]) + fabs(bacm[2][3] - bacm[4][3])) / 3;
    v1v2c = (v1v2c / 3); // - 0.5) * 2;
    // R1-R4 and R2-R3 horizontal simmetry
    float r1r4h = (fabs(bacm[0][0] - bacm[6][0]) + fabs(bacm[0][1] - bacm[6][1]) + fabs(bacm[0][2] - bacm[6][2]) +
                   fabs(bacm[1][0] - bacm[5][0]) + fabs(bacm[1][1] - bacm[5][1]) + fabs(bacm[1][2] - bacm[5][2]) +
                   fabs(bacm[2][0] - bacm[4][0]) + fabs(bacm[2][1] - bacm[4][1]) + fabs(bacm[2][2] - bacm[4][2])) / 9;
    r1r4h = (r1r4h / 9); // - 0.5) * 2;
    float r2r3h = (fabs(bacm[0][4] - bacm[6][4]) + fabs(bacm[0][5] - bacm[6][5]) + fabs(bacm[0][6] - bacm[6][6]) +
                   fabs(bacm[1][4] - bacm[5][4]) + fabs(bacm[1][5] - bacm[5][5]) + fabs(bacm[1][6] - bacm[5][6]) +
                   fabs(bacm[2][4] - bacm[4][4]) + fabs(bacm[2][5] - bacm[4][5]) + fabs(bacm[2][6] - bacm[4][6])) / 9;
    r2r3h = (r2r3h / 9); // - 0.5) * 2;
    // R1-R2 and R4-R3 vertical simmetry
    float r1r2v = (fabs(bacm[0][0] - bacm[0][6]) + fabs(bacm[0][1] - bacm[0][5]) + fabs(bacm[0][2] - bacm[0][4]) +
                   fabs(bacm[1][0] - bacm[1][6]) + fabs(bacm[1][1] - bacm[1][5]) + fabs(bacm[1][2] - bacm[1][4]) +
                   fabs(bacm[2][0] - bacm[2][6]) + fabs(bacm[2][1] - bacm[2][5]) + fabs(bacm[2][2] - bacm[2][4])) / 9;
    r1r2v = (r1r2v / 9); // - 0.5) * 2;
    float r4r3v = (fabs(bacm[4][0] - bacm[4][6]) + fabs(bacm[4][1] - bacm[4][5]) + fabs(bacm[4][2] - bacm[4][4]) +
                   fabs(bacm[5][0] - bacm[5][6]) + fabs(bacm[5][1] - bacm[5][5]) + fabs(bacm[5][2] - bacm[5][4]) +
                   fabs(bacm[6][0] - bacm[6][6]) + fabs(bacm[6][1] - bacm[6][5]) + fabs(bacm[6][2] - bacm[6][4])) / 9;
    r4r3v = (r4r3v / 9); // - 0.5) * 2;
    // R1-R3 and R2-R4 central simmetry
    float r1r3c = (fabs(bacm[0][0] - bacm[6][6]) + fabs(bacm[0][1] - bacm[6][5]) + fabs(bacm[0][2] - bacm[6][4]) +
                   fabs(bacm[1][0] - bacm[5][6]) + fabs(bacm[1][1] - bacm[5][5]) + fabs(bacm[1][2] - bacm[5][4]) +
                   fabs(bacm[2][0] - bacm[4][6]) + fabs(bacm[2][1] - bacm[4][5]) + fabs(bacm[2][2] - bacm[4][4])) / 9;
    r1r3c = (r1r3c / 9); // - 0.5) * 2;
    float r2r4c = (fabs(bacm[0][4] - bacm[6][2]) + fabs(bacm[0][5] - bacm[6][1]) + fabs(bacm[0][6] - bacm[6][0]) +
                   fabs(bacm[1][4] - bacm[5][2]) + fabs(bacm[1][5] - bacm[5][1]) + fabs(bacm[1][6] - bacm[5][0]) +
                   fabs(bacm[2][4] - bacm[4][2]) + fabs(bacm[2][5] - bacm[4][1]) + fabs(bacm[2][6] - bacm[4][0])) / 9;
    r2r4c = (r2r4c / 9); // - 0.5) * 2;

    // Percentage of C against H, V, R1, R2, R3, R4
    float hh_avg = (bacm[3][0] + bacm[3][1] + bacm[3][2] + bacm[3][4] + bacm[3][5] + bacm[3][6]); // / 6;
    float vv_avg = (bacm[0][3] + bacm[1][3] + bacm[2][3] + bacm[4][3] + bacm[5][3] + bacm[6][3]); // / 6;
    float r1_avg = 0;
    for (int m = 0; m <= 2; ++m)
        for (int n = 0; n <= 2; ++n)
            r1_avg += bacm[m][n];
    // r1_avg = r1_avg / 9 + numeric_limits<float>::epsilon();
    float r2_avg = 0;
    for (int m = 0; m <= 2; ++m)
        for (int n = 4; n <= 6; ++n)
            r2_avg += bacm[m][n];
    // r2_avg = r2_avg / 9 + numeric_limits<float>::epsilon();
    float r3_avg = 0;
    for (int m = 4; m <= 6; ++m)
        for (int n = 4; n <= 6; ++n)
            r3_avg += bacm[m][n];
    // r3_avg = r3_avg / 9 + numeric_limits<float>::epsilon();
    float r4_avg = 0;
    for (int m = 4; m <= 6; ++m)
        for (int n = 0; n <= 2; ++n)
            r4_avg += bacm[m][n];
    // r4_avg = r4_avg / 9 + numeric_limits<float>::epsilon();

    float c = bacm[3][3];
    float cr1 = r1_avg != 0 ? (c / r1_avg) / 16 : 0;
    float cr2 = r2_avg != 0 ? (c / r2_avg) / 16 : 0;
    float cr3 = r3_avg != 0 ? (c / r3_avg) / 16 : 0;
    float cr4 = r4_avg != 0 ? (c / r4_avg) / 16 : 0;
    float chh = hh_avg != 0 ? (c / hh_avg) / 16 : 0;
    float cvv = vv_avg != 0 ? (c / vv_avg) / 16 : 0;

    features.create(1, NUM_FEATURES, CV_32F);
    features.at<float>(0, 0) = h1h2c;
    features.at<float>(0, 1) = v1v2c;
    features.at<float>(0, 2) = r1r4h;
    features.at<float>(0, 3) = r2r3h;
    features.at<float>(0, 4) = r1r2v;
    features.at<float>(0, 5) = r4r3v;
    features.at<float>(0, 6) = r1r3c;
    features.at<float>(0, 7) = r2r4c;
    features.at<float>(0, 8) = cr1;
    features.at<float>(0, 9) = cr2;
    features.at<float>(0, 10) = cr3;
    features.at<float>(0, 11) = cr4;
    features.at<float>(0, 12) = chh;
    features.at<float>(0, 13) = cvv;

    // features.clear();
    // features.push_back(h1h2c);
    // features.push_back(v1v2c);
    // features.push_back(r1r4h);
    // features.push_back(r2r3h);
    // features.push_back(r1r2v);
    // features.push_back(r4r3v);
    // features.push_back(r1r3c);
    // features.push_back(r2r4c);
    // features.push_back(cr1);
    // features.push_back(cr2);
    // features.push_back(cr3);
    // features.push_back(cr4);
    // features.push_back(chh);
    // features.push_back(cvv);

    if (DEBUG)
    {
        printf("\n\nListContourPlot[{");
        for (int m = 0; m < BLOCK_SIZE; ++m)
        {
            for (int n = 0; n < BLOCK_SIZE; ++n)
            {
                printf("{%d,%d,%f}", m + 1, n + 1, bacm[m][n]);
                if (m != BLOCK_SIZE - 1 || n != BLOCK_SIZE - 1)
                    printf(",");
            }
        }
        printf("}]\n");
        printf("ListPlot[{");
        for (int k = 0; k < NUM_FEATURES; ++k)
        {
            printf("%f", features.at<float>(0, k));
            if (k < NUM_FEATURES - 1)
                printf(",");
        }
        printf("}, Filling->Axis, PlotRange->All]\n");
    }
}

void writeFeatures(FILE* file, int label, Mat features)
{
    string line = to_string(label) + " ";
    for (int i = 0; i < features.cols; ++i)
        line += to_string(i + 1) + ":" + to_string(features.at<float>(0, i)) + " ";
    line += "\n";
    fprintf(file, "%s", line.c_str());
}

Mat randomCrop(const Mat input)
{
    int crop_row = rand() % 8;
    int crop_col = rand() % 8;
    Mat cropped(input, Rect(crop_col, crop_row, input.cols - crop_col, input.rows - crop_row));
    return (cropped.clone());
}

void createTraining(bool crop = false)
{
    const string PATH("/home/vision/Downloads/ucid/");
    const string TRAIN_FILE1(crop ? "svm/train_c.txt" : "svm/train.txt");
    const string TRAIN_FILE2(crop ? "svm/train_c.yml" : "svm/train.yml");
    DIR* dir = opendir(PATH.c_str());
    struct dirent *ent;
    vector<string> filenames;
    while ((ent = readdir (dir)) != NULL)
        filenames.push_back(string(ent->d_name));
    closedir(dir);
    int total_files = (int)filenames.size();

    const int START_QUAL = 50;
    const int END_QUAL = 90;
    const int Q1_STEP = 5;
    const int Q2_STEP = 5;
    const int NUM_FEATURES = 180;
    const int NUM_Q1 = (END_QUAL - START_QUAL) / Q1_STEP + 1;
    const int NUM_Q2 = (END_QUAL - START_QUAL) / Q2_STEP + 1;
    const int REPEAT_Q1 = (int)(total_files / (float)NUM_Q1);
    const int REPEAT_Q2 = (int)(total_files / (float)(NUM_Q1 * NUM_Q2));
    const int TRAIN_SIZE = REPEAT_Q1 * NUM_Q1 + REPEAT_Q2 * NUM_Q1 * NUM_Q2;

    Mat train_in(TRAIN_SIZE, NUM_FEATURES, CV_32F, Scalar(0));
    Mat train_out(TRAIN_SIZE, 1, CV_32S, Scalar(-1));

    uint train_index = 0;
    uint file_index = 0;
    FILE* train_file = fopen(TRAIN_FILE1.c_str(), "w");
    FileStorage fs(TRAIN_FILE2, FileStorage::WRITE);
    for (int r = 0; r < REPEAT_Q1; ++r)
    {
        for (int q1 = START_QUAL; q1 <= END_QUAL; q1 += Q1_STEP)
        {
            string filename = filenames.at(file_index);
            Mat image = imread(PATH + filename, CV_LOAD_IMAGE_COLOR);
            if (!image.empty())
            {
                printf("{%d/%d} %d/%d: '%s' (%d%%)\n", r + 1, REPEAT_Q1, file_index + 1, total_files, filename.c_str(), q1);
                Mat q1_comp = jpegCompress(image, q1);
                Mat features;
                firstDigitFeatures(q1_comp, features);
                // bacmFeatures(q1_comp, features);
                writeFeatures(train_file, 0, features);
                for (int k = 0; k < NUM_FEATURES; ++k)
                    train_in.at<float>(train_index, k) = features.at<float>(0, k);
                train_out.at<int>(train_index, 0) = 0;
                ++train_index;
            }
            if (++file_index == filenames.size())
                file_index = 0;
        }
    }

    file_index = 0;
    for (int r = 0; r < REPEAT_Q2; ++r)
    {
        for (int q1 = START_QUAL; q1 <= END_QUAL; q1 += Q1_STEP)
        {
            for (int q2 = START_QUAL; q2 <= END_QUAL; q2 += Q2_STEP)
            {
                if (q1 == q2)
                    continue;
                string filename = filenames.at(file_index);
                Mat image = imread(PATH + filename, CV_LOAD_IMAGE_COLOR);
                if (!image.empty())
                {
                    printf("{%d/%d} %d/%d: '%s' (%d%% -> %d%%)\n", r + 1, REPEAT_Q2, file_index + 1, total_files, filename.c_str(), q1, q2);
                    Mat q1_comp = jpegCompress(image, q1);
                    Mat q1_cropped;
                    if (crop)
                        q1_cropped = randomCrop(q1_comp);
                    else
                        q1_cropped = q1_comp;
                    Mat q2_comp = jpegCompress(q1_cropped, q2);
                    Mat features;
                    firstDigitFeatures(q2_comp, features);
                    // bacmFeatures(q2_comp, features);
                    writeFeatures(train_file, 1, features);
                    for (int k = 0; k < NUM_FEATURES; ++k)
                        train_in.at<float>(train_index, k) = features.at<float>(0, k);
                    train_out.at<int>(train_index, 0) = 1;
                    ++train_index;
                }
                if (++file_index == filenames.size())
                    file_index = 0;
            }
        }
    }

    fs << "train_in" << train_in.rowRange(0, train_index);
    fs << "train_out" << train_out.rowRange(0, train_index);
    fclose(train_file);
}



void trainSVM(bool crop = false)
{
    const string TRAIN_FILE(crop ? "svm/train_c.yml" : "svm/train.yml");
    const string SVM_FILE(crop ? "svm/model_c.yml" : "svm/model.yml");

    FileStorage fs(TRAIN_FILE, FileStorage::READ);
    Mat train_in, train_out;
    fs["train_in"] >> train_in;
    fs["train_out"] >> train_out;

    CvSVMParams params;
    params.svm_type    = CvSVM::C_SVC;
    params.kernel_type = CvSVM::RBF;
    params.gamma       = 0.0625;
    params.C           = 64;
    params.term_crit   = cvTermCriteria(CV_TERMCRIT_EPS, 1500000, 1e-3);
    CvSVM svm;
    // svm.train_auto(train_in, train_out, Mat(), Mat(), params,
    //                 5, // Validation folds
    //                 CvParamGrid(1, 1024, 2), // C
    //                 CvParamGrid(0.005, 0.5, 0), // G
    //                 CvSVM::get_default_grid(CvSVM::P),
    //                 CvSVM::get_default_grid(CvSVM::NU),
    //                 CvSVM::get_default_grid(CvSVM::COEF),
    //                 CvSVM::get_default_grid(CvSVM::DEGREE),
    //                 true);
    svm.train(train_in, train_out, Mat(), Mat(), params);
    printf("C = %f, G = %f, N = %d\n", svm.get_params().C, svm.get_params().gamma, svm.get_support_vector_count());
    svm.save(SVM_FILE.c_str());
}

int dctHistogram(const Mat image, float &probability)
{
    const int   SKIP_DC   = 1;
    const int   NUM_MODES = 9;
    const short MIN_COEFF = -128;
    const short MAX_COEFF = +128;
    const short NUM_BINS  = MAX_COEFF - MIN_COEFF + 1;

    Mat gray;
    cvtColor(image, gray, CV_BGR2GRAY);
    gray.convertTo(gray, CV_32F);
    Mat coeff_hist(NUM_MODES, NUM_BINS, CV_32S, Scalar(0));
    vector<int> coeff_values[NUM_MODES];
    for (int i = 0; i < gray.rows - DCT_SIZE; i += DCT_SIZE)
    {
        for (int j = 0; j < gray.cols - DCT_SIZE; j += DCT_SIZE)
        {
            Mat block(gray, Rect(j, i, DCT_SIZE, DCT_SIZE));
            block -= DCT_OFFSET;
            Mat coeffs(DCT_SIZE, DCT_SIZE, CV_32F);
            dct(block, coeffs);
            coeffs.convertTo(coeffs, CV_16S);
            for (int m = 0; m < NUM_MODES; ++m)
            {
                int i0 = ZIG_ZAG.at<int>(m + SKIP_DC, 0);
                int j0 = ZIG_ZAG.at<int>(m + SKIP_DC, 1);
                short coeff = coeffs.at<short>(i0, j0);
                coeff_values[m].push_back(coeff);
                if (coeff >= MIN_COEFF && coeff <= MAX_COEFF)
                    ++coeff_hist.at<int>(m, coeff - MIN_COEFF);
            }
        }
    }

    // for (int m = 0; m < NUM_MODES; ++m)
    // {
    //     float coeff_sum = 0;
    //     for (int b = 0; b < NUM_BINS; ++b)
    //         coeff_sum += coeff_hist.at<float>(m, b);
    //     for (int b = 0; b < NUM_BINS; ++b)
    //         coeff_hist.at<float>(m, b) /= coeff_sum;
    // }

    // for (uint m = 0; m < NUM_MODES; ++m)
    // {
    //     printf("hist%d={", m);
    //     for (uint b = 0; b < NUM_BINS; ++b)
    //         printf("%d%s", coeff_hist.at<int>(m, b), b < NUM_BINS - 1? "," : "");
    //     printf("};\n");
    // }

    for (uint m = 0; m < NUM_MODES; ++m)
    {
        if (m != 3)
            continue;
        printf("hist%d={", m);
        for (uint c = 0; c < coeff_values[m].size(); ++c)
            printf("%d%s", coeff_values[m].at(c), c < coeff_values[m].size() - 1? "," : "");
        printf("};\n");
    }

    return 0;
}

void analyzeProbability()
{
    const string SVM_FILE("svm/model.yml");
    const string PATH("/home/vision/Downloads/ucid/");
    DIR* dir = opendir(PATH.c_str());
    struct dirent *ent;
    vector<string> filenames;
    while ((ent = readdir (dir)) != NULL)
        filenames.push_back(string(ent->d_name));
    closedir(dir);
    int total_files = (int)filenames.size();

    CvSVM svm;
    svm.load(SVM_FILE.c_str());
    vector<float> probs;

    const int START_QUAL = 50;
    const int END_QUAL = 90;
    const int Q1_STEP = 5;
    const int Q2_STEP = 5;
    const int NUM_Q1 = (END_QUAL - START_QUAL) / Q1_STEP + 1;
    const int NUM_Q2 = (END_QUAL - START_QUAL) / Q2_STEP + 1;
    const int REPEAT_Q1 = (int)(total_files / (float)NUM_Q1);
    const int REPEAT_Q2 = (int)(total_files / (float)(NUM_Q1 * NUM_Q2));
    const int FILE_LIMIT = 700;

    bool stop = false;
    uint file_index = 0;
    for (int r = 0; r < REPEAT_Q1; ++r)
    {
        for (int q1 = START_QUAL; q1 <= END_QUAL; q1 += Q1_STEP)
        {
            string filename = filenames.at(file_index);
            Mat image = imread(PATH + filename, CV_LOAD_IMAGE_COLOR);
            if (!image.empty())
            {
                printf("{%d/%d} %d/%d: '%s' (%d%%)\n", r + 1, REPEAT_Q1, file_index + 1, total_files, filename.c_str(), q1);
                Mat q1_comp = jpegCompress(image, q1);
                Mat features;
                firstDigitFeatures(q1_comp, features);
                probs.push_back(svm.predict(features, true));
            }
            if (++file_index == FILE_LIMIT)
                stop = true;
            if (stop)
                break;
        }
        if (stop)
            break;
    }

    file_index = 0;
    stop = false;
    for (int r = 0; r < REPEAT_Q2; ++r)
    {
        for (int q1 = START_QUAL; q1 <= END_QUAL; q1 += Q1_STEP)
        {
            for (int q2 = START_QUAL; q2 <= END_QUAL; q2 += Q2_STEP)
            {
                if (q1 == q2)
                    continue;
                string filename = filenames.at(file_index);
                Mat image = imread(PATH + filename, CV_LOAD_IMAGE_COLOR);
                if (!image.empty())
                {
                    printf("{%d/%d} %d/%d: '%s' (%d%% -> %d%%)\n", r + 1, REPEAT_Q2, file_index + 1, total_files, filename.c_str(), q1, q2);
                    Mat q1_comp = jpegCompress(image, q1);
                    Mat q2_comp = jpegCompress(q1_comp, q2);
                    Mat features;
                    firstDigitFeatures(q2_comp, features);
                    probs.push_back(svm.predict(features, true));
                }
                if (++file_index == FILE_LIMIT)
                stop = true;
                if (stop)
                    break;
            }
            if (stop)
                break;
        }
        if (stop)
            break;
    }

    printf("data={");
    for (uint i = 0; i < probs.size(); ++i)
        printf("%f,", probs.at(i));
    printf("};\n");
}

void testLossless()
{
    const string PATH("/home/vision/Downloads/ucid/");
    DIR* dir = opendir(PATH.c_str());
    struct dirent *ent;
    vector<string> filenames;
    while ((ent = readdir (dir)) != NULL)
        filenames.push_back(string(ent->d_name));
    closedir(dir);
    filenames.erase(filenames.begin(), filenames.begin() + 2);
    int total_files = (int)filenames.size();

    const int START_QUAL = 10;
    const int END_QUAL = 90;
    const int Q1_STEP = 1;
    const int NUM_Q1 = (END_QUAL - START_QUAL) / Q1_STEP + 1;
    const int REPEAT_Q1 = (int)(total_files / (float)NUM_Q1);
    const int FILE_LIMIT = 1300;

    FileStorage fs("peaks.yml", FileStorage::WRITE);
    Mat peaks1(1, FILE_LIMIT, CV_32F, Scalar(0));
    Mat peaks2(1, FILE_LIMIT, CV_32F, Scalar(0));

    bool stop = false;
    uint file_index = 0;
    uint peak_index = 0;
    for (int r = 0; r < REPEAT_Q1; ++r)
    {
        for (int q1 = START_QUAL; q1 <= END_QUAL; q1 += Q1_STEP)
        {
            string filename = filenames.at(file_index);
            Mat image = imread(PATH + filename, CV_LOAD_IMAGE_COLOR);
            if (!image.empty())
            {
                printf("{%d/%d} %d/%d: '%s' (%d%%)\n", r + 1, REPEAT_Q1, file_index + 1, total_files, filename.c_str(), q1);
                int l1, l2, n, q1, q2;
                float p1, p2;
                Mat plot;
                compressionGhosts(image, n, l1, p1, l2, p2, q1, q2, plot);
                peaks1.at<float>(0, peak_index) = p1;
                Mat comp = jpegCompress(image, q1);
                compressionGhosts(comp, n, l1, p1, l2, p2, q1, q2, plot);
                peaks2.at<float>(0, peak_index) = p1;
                ++peak_index;
            }
            if (++file_index == FILE_LIMIT)
            {
                stop = true;
                break;
            }
        }
        if (stop)
            break;
    }

    fs << "peaks1" << peaks1.colRange(0, peak_index);
    fs << "peaks2" << peaks2.colRange(0, peak_index);

    printf("peaks1={");
    for (uint i = 0; i < peak_index; ++i)
        printf("%f,", peaks1.at<float>(0, i));
    printf("};\n");
    printf("peaks2={");
    for (uint i = 0; i < peak_index; ++i)
        printf("%f,", peaks2.at<float>(0, i));
    printf("};\n");
    printf("Histogram[{peaks1,peaks2}]\n");
}

void testLossless2()
{
    FileStorage fs("peaks.yml", FileStorage::READ);
    Mat peaks1, peaks2;
    fs["peaks1"] >> peaks1;
    fs["peaks2"] >> peaks2;

    float avg = 0;
    for (int i = 0; i < peaks1.cols; ++i)
        avg += peaks1.at<float>(0, i) + peaks2.at<float>(0, i);
    avg /= peaks1.cols * 2;
    printf("T0 = %f\n", avg);

    float dist = 0;
    for (int i = 0; i < peaks1.cols; ++i)
        dist += (peaks1.at<float>(0, i) + peaks2.at<float>(0, i)) / 2;
    dist /= peaks1.cols;
    printf("T1 = %f\n", dist);

    const float THR_EPS = 0.00001;
    float thr_prev = 1.0;
    float thr_curr = 0.5;
    float a = 0;
    float b = 1;
    int iter = 0;
    while (fabs(thr_curr - thr_prev) > THR_EPS)
    {
        float dist1 = 0;
        float dist2 = 0;
        for (int i = 0; i < peaks1.cols; ++i)
        {
            dist1 += fabs(peaks1.at<float>(0, i) - thr_curr);
            dist2 += fabs(peaks2.at<float>(0, i) - thr_curr);
        }
        if (dist1 > dist2)
            b = thr_curr;
        else
            a = thr_curr;
        thr_prev = thr_curr;
        thr_curr = (a + b) / 2;
        printf("%d) T2 = %f\n", ++iter, thr_curr);
    }
}

// ----------------------------------------------------------------------------------- //

void firstDigitPDF(const Mat input, int quality, vector<float> &digits_pdf)
{
    const int BLOCK_SIZE = 8;
    const int MAX_DIGITS = 9;
    const bool DEBUG = false;

    Mat gray;
    cvtColor(input, gray, CV_BGR2GRAY);
    gray.convertTo(gray, CV_32F);
    Mat luma_qt, chroma_qt, table;
    getQTs(quality, luma_qt, chroma_qt);
    // printMat("LUMA_QT", luma_qt);
    luma_qt.convertTo(table, CV_32F);
    digits_pdf = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<float> pdf_min = {0.248635, 0.037269, 0.020296, 0.010299, 0.005706, 0.003754, 0.002689, 0.000894, 0.003698};
    vector<float> pdf_max = {0.847281, 0.270014, 0.200408, 0.169191, 0.11372, 0.183909, 0.125889, 0.108814, 0.150106};

    for (int i = 0; i < input.rows - BLOCK_SIZE; i += BLOCK_SIZE)
    {
        for (int j = 0; j < input.cols - BLOCK_SIZE; j += BLOCK_SIZE)
        {
            Mat block(gray, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            // printMat("BLOCK1", block);
            block -= 128;
            // printMat("BLOCK2", block);
            Mat coeffs(BLOCK_SIZE, BLOCK_SIZE, CV_32F);
            dct(block, coeffs);
            // printMat("COEFFS", coeffs);
            coeffs /= table;
            // printMat("COEFFS_Q1", coeffs);
            coeffs.convertTo(coeffs, CV_16S);
            // printMat("COEFFS_Q2", coeffs);
            for (int h = 0; h < BLOCK_SIZE; ++h)
            {
                for (int k = 0; k < BLOCK_SIZE; ++k)
                {
                    if (h == 0 && k == 0)
                        continue;
                    int digit = getFirstDigit(coeffs.at<short>(h, k));
                    if (digit > 0)
                        ++digits_pdf[digit - 1];
                }
            }
        }
    }

    float digits_sum = 0;
    for (int i = 0; i < MAX_DIGITS; ++i)
        digits_sum += digits_pdf.at(i);
    for (int i = 0; i < MAX_DIGITS; ++i)
    {
        digits_pdf[i] /= digits_sum;
        // digits_pdf[i] -= pdf_min.at(i);
        // digits_pdf[i] /= (pdf_max.at(i) - pdf_min.at(i));
        digits_pdf[i] *= 2;
        digits_pdf[i] -= 1;
    }

    if (DEBUG)
    {
        printf("{");
        for (int i = 0; i < MAX_DIGITS; ++i)
            printf("%f,", digits_pdf.at(i));
        printf("}");
    }
}

void doubleCompression4(const Mat input, float &probability)
{
    const string FILENAME("libsvm/tools/test.txt");
    Mat features;
    // printf("\ndata = ");
    firstDigitFeatures(input, features);
    // printf(";\n");
    FILE* test_file = fopen(FILENAME.c_str(), "w");
    writeFeatures(test_file, 1, features);
    fclose(test_file);

    // const int MIN_QUALITY = 50;
    // const int MAX_QUALITY = 85;
    // const int QUALITY_STEP = 5;

    // printf("\ndata = {");
    // for (int q = MIN_QUALITY; q <= MAX_QUALITY; q += QUALITY_STEP)
    // {
    //     vector<float> pdf;
    //     firstDigitPDF(input, q, pdf);
    // }
    // printf("};\n");
    // printf("ListLinePlot[data, PlotLegends->{");
    // for (int q = MIN_QUALITY; q <= MAX_QUALITY; q += QUALITY_STEP)
    //     printf("\"%d\",", q);
    // printf("}, PlotRange->All, PlotStyle->Thick]\n");
}

void doubleCompression9(const Mat input, float &probability)
{
    const int TOTAL_FILES = 1375;
    const int START_QUAL = 50;
    const int END_QUAL = 100;
    const int QUAL_STEP = 10;
    const string PATH("/home/vision/Downloads/ucid/");

    DIR* dir = opendir(PATH.c_str());
    struct dirent *ent;
    vector<string> filenames;
    while ((ent = readdir (dir)) != NULL)
        filenames.push_back(string(ent->d_name));
    closedir(dir);

    vector< vector<float> > pdfs;
    for (int q = START_QUAL; q <= END_QUAL; q += QUAL_STEP)
    {
        vector<float> pdf_qual = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        for (uint f = 0; f <= filenames.size() - 2; ++f)
        {
            string filename = filenames.at(f);
            Mat image = imread(PATH + filename, CV_LOAD_IMAGE_COLOR);
            if (!image.empty())
            {
                Mat comp = jpegCompress(image, q);
                vector<float> pdf;
                firstDigitPDF(comp, q, pdf);
                for (uint i = 0; i < pdf.size(); ++i)
                    pdf_qual[i] += pdf.at(i);
                printf("%d/%d: '%s' (%d%%)\n", f + 1, TOTAL_FILES, filename.c_str(), q);
            }
        }
        float pdf_sum = 0;
        for (uint i = 0; i < pdf_qual.size(); ++i)
            pdf_sum += pdf_qual.at(i);
        for (uint i = 0; i < pdf_qual.size(); ++i)
            pdf_qual[i] /= pdf_sum;
        pdfs.push_back(pdf_qual);
    }

    for (uint i = 0; i < pdfs.size(); ++i)
    {
        printf("%d = {", i);
        for (uint j = 0; j < pdfs.at(i).size(); ++j)
            printf("%f,", pdfs.at(i).at(j));
        printf("}\n");
    }
}

void doubleCompression10(const Mat input, float &probability)
{
    const string FILENAME("test.txt");
    vector<float> features;
    firstDigitPDF(input, 100, features);
    FILE* test_file = fopen(FILENAME.c_str(), "w");
    // writeFeatures(test_file, 0, features);
    fclose(test_file);
}

// void doubleCompression4(const Mat input, float &probability)
// {
//     const int TOTAL_FILES = 1375;
//     const int START_QUAL = 50;
//     const int END_QUAL = 95;
//     const int Q2_STEP = 1;
//     const int NUM_FEATURES = 14;
//     const string PATH("/home/vision/Downloads/ucid/");
//     const string FILENAME("train1.txt");

//     FILE* train_file = fopen(FILENAME.c_str(), "w");
//     DIR* dir = opendir(PATH.c_str());
//     struct dirent *ent;
//     int file_count = 0;

//     vector< vector<int> > factors;
//     for (int q1 = START_QUAL; q1 <= END_QUAL; ++q1)
//     {
//         for (int q2 = START_QUAL; q2 <= END_QUAL; q2 += Q2_STEP)
//         {
//             vector<int> q1q2;
//             q1q2.push_back(q1);
//             q1q2.push_back(q2);
//             factors.push_back(q1q2);
//         }
//     }
//     int num_q1 = END_QUAL - START_QUAL + 1;
//     int num_q2 = (END_QUAL - START_QUAL) / Q2_STEP + 1;
//     // printf("N1 = %d, N2 = %d\n", num_q1, num_q2);
//     Mat train_data(num_q1 + num_q1 * num_q2, NUM_FEATURES, CV_32F);
//     Mat labels(train_data.rows, 1, CV_32F);
//     // printf("ROWS = %d, COLS = %d", train_data.rows, train_data.cols);
//     int train_count = 0;

//     printf("\n");
//     while ((ent = readdir (dir)) != NULL)
//     {
//         if (file_count == (int)factors.size())
//             break;
//         Mat image = imread(PATH + string(ent->d_name), CV_LOAD_IMAGE_COLOR);
//         if (image.empty())
//             continue;

//         int q1 = factors.at(file_count).at(0);
//         int q2 = factors.at(file_count).at(1);
//         printf("%d/%d: '%s' ", file_count + 1, TOTAL_FILES, ent->d_name);
//         vector<float> features;

//         printf("(%d%%) >> [", q1);
//         Mat q1_comp = jpegCompress(image, q1);
//         if (file_count == 0 || q1 != factors.at(file_count - 1).at(0))
//         {
//             bacmFeatures(q1_comp, features);
//             writeFeatures(train_file, 0, features);
//             for (int i = 0; i < NUM_FEATURES; ++i)
//                 train_data.at<float>(train_count, i) = features.at(i);
//             labels.at<float>(train_count, 0) = 0;
//             ++train_count;
//         }
//         printf("] ");

//         Mat q1_cropped = randomCrop(q1_comp);
//         Mat q2_comp = jpegCompress(q1_cropped, q2);
//         printf("(%d%%) >> [", q2);
//         bacmFeatures(q2_comp, features);
//         printf("]\n");
//         writeFeatures(train_file, 1, features);
//         for (int i = 0; i < NUM_FEATURES; ++i)
//             train_data.at<float>(train_count, i) = features.at(i);
//         labels.at<float>(train_count, 0) = 1;
//         ++train_count;
//         ++file_count;
//     }
//     closedir(dir);
//     fclose(train_file);

//     FileStorage fs("train1.xml", FileStorage::WRITE);
//     fs << "train" << train_data;
//     fs << "labels" << labels;
// }

// void doubleCompression5(const Mat input, float &probability)
// {
//     vector<float> features;
//     bacmFeatures(input, features);
// }

// void doubleCompression6(const Mat input, float &probability)
// {
//     FILE* test_file = fopen("test.txt", "w");
//     FileStorage fs("train2.xml", FileStorage::READ);
//     Mat train_data, labels;
//     fs["train"] >> train_data;
//     fs["labels"] >> labels;

//     CvSVMParams params;
//     params.svm_type    = CvSVM::C_SVC;
//     params.kernel_type = CvSVM::RBF;
//     params.gamma       = 0.015625;
//     params.C           = 0.03125;
//     params.term_crit   = cvTermCriteria(CV_TERMCRIT_EPS, 100, 1e-6);
//     CvSVM svm;
//     svm.train(train_data, labels, Mat(), Mat(), params);

//     vector<float> features;
//     bacmFeatures(input, features);
//     Mat sample(1, 14, CV_32F);
//     for (int i = 0; i < 14; ++i)
//         sample.at<float>(0, i) = features.at(i);
//     printf("P = %f\n", svm.predict(sample, false));
//     writeFeatures(test_file, -1, features);

//     fclose(test_file);
// }

int computeBlock(const Mat mat, int row, int col)
{
    if (row < 0 || row > mat.rows - 2 || col < 0 || col > mat.cols - 2)
        return 0;
    int a = mat.at<uchar>(row, col);
    int b = mat.at<uchar>(row, col + 1);
    int c = mat.at<uchar>(row + 1, col);
    int d = mat.at<uchar>(row + 1, col + 1);
    return abs(a - b - c + d);
}

void gridDetection(const Mat input, int &x, int &y)
{
    const int BLOCK_SIZE = 8;

    Mat gray;
    cvtColor(input, gray, CV_BGR2GRAY);
    Mat err(BLOCK_SIZE, BLOCK_SIZE, CV_32F, Scalar(0));
    for (int i0 = 0; i0 < BLOCK_SIZE; ++i0)
    {
        int sum = 0;
        for (int i = i0 + BLOCK_SIZE - 1; i < gray.rows - BLOCK_SIZE; i += BLOCK_SIZE)
            for (int j = 0; j < gray.cols; ++j)
                sum += computeBlock(gray, i, j);
        err.at<float>(i0, 0) = sum;
    }
    printMat("ERR", err, true);
}

void kurtosisAnalysis(const Mat input, Mat &output)
{
    return;
    const int BLOCK_SIZE = 16;
    const int BLOCK_ELEM = BLOCK_SIZE * BLOCK_SIZE;

    Mat gray;
    cvtColor(input, gray, CV_BGR2GRAY);
    output.create(input.rows, input.cols, CV_32F);
    output.setTo(0);
    for (int i = BLOCK_SIZE; i < input.rows - BLOCK_SIZE; ++i)
    {
        for (int j = BLOCK_SIZE; j < input.cols - BLOCK_SIZE; ++j)
        {
            Mat block(gray, Rect(j - BLOCK_SIZE, i - BLOCK_SIZE, BLOCK_SIZE, BLOCK_SIZE));
            float avg = mean(block)[0];
            float den = 0;
            float num = 0;
            for (int h = 0; h < BLOCK_SIZE; ++h)
            {
                for (int k = 0; k < BLOCK_SIZE; ++k)
                {
                    float dev = (int)block.at<uchar>(h, k) - avg;
                    num += pow(dev, 4);
                    den += pow(dev, 2);
                }
            }
            float kurt = (num*BLOCK_ELEM*(BLOCK_ELEM+1)*(BLOCK_ELEM-1)/(den*den*(BLOCK_ELEM-2)*(BLOCK_ELEM-3)))-(3*(BLOCK_ELEM-1)*(BLOCK_ELEM-1)/((BLOCK_ELEM-2)*(BLOCK_ELEM-3)));
            output.at<float>(i, j) = kurt;
            // Mat block2(output, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            // block2 = kurt;
        }
    }
    normalize(output, output, 0, 255, NORM_MINMAX);
    output.convertTo(output, CV_8U);
    equalizeHist(output, output);
}

double findThreshold(const Mat m1, const Mat m2)
{
    const double STEP = 0.001;

    double min1, max1, min2, max2;
    minMaxLoc(m1, &min1, &max1);
    minMaxLoc(m2, &min2, &max2);
    double min0 = min(min1, min2);
    double max0 = max(max1, max2);
    int num_steps = (int)((max0 - min0) / STEP) + 1;

    double min_err = numeric_limits<double>::max();
    double best_t = -1;
    for (int i = 0; i < num_steps; ++i)
    {
        double t = min0 + i * STEP;
        int err = 0;
        for (int i = 0; i < m1.cols; ++i)
            if (m1.at<float>(0, i) > t)
                ++err;
        for (int i = 0; i < m2.cols; ++i)
            if (m2.at<float>(0, i) < t)
                ++err;
        if (err < min_err)
        {
            min_err = err;
            best_t = t;
        }
    }
    return best_t;
}

void testContrast()
{
    const uchar MIN_LOW  = 50;
    const uchar MAX_LOW  = 120;
    const uchar MIN_HIGH = 130;
    const uchar MAX_HIGH = 200;
    const string PATH("/home/bartoli/Downloads/ucid/");

    DIR* dir = opendir(PATH.c_str());
    struct dirent *ent;
    vector<string> filenames;
    while ((ent = readdir (dir)) != NULL)
        filenames.push_back(string(ent->d_name));
    closedir(dir);
    filenames.erase(filenames.begin(), filenames.begin() + 2);
    int num_files = (int)filenames.size();

    Mat errors1(1, num_files, CV_32F);
    Mat errors2(1, num_files, CV_32F);
    for (int i = 0; i < num_files; ++i)
    {
        string filename = filenames.at(i);
        Mat image = imread(PATH + filename, CV_LOAD_IMAGE_COLOR);
        double err1 = computeChanSim(image);
        errors1.at<float>(0, i) = err1;

        uchar low = (uchar)(rand() % (MAX_LOW - MIN_LOW) + MIN_LOW);
        uchar high = (uchar)(rand() % (MAX_HIGH - MIN_HIGH) + MIN_HIGH);
        Mat lut = buildContrastLUT(low, high);
        Mat image2;
        LUT(image, lut, image2);
        double err2 = computeChanSim(image2);
        errors2.at<float>(0, i) = err2;

        printf("[%d/%d, '%s',{%d,%d}]:\tE1 = %f, E2 = %f\n", i + 1, num_files, filename.c_str(), low, high, err1, err2);
    }

    printf("e1=[");
    for (int i = 0; i < errors1.cols; ++i)
        printf("%f,", errors1.at<float>(0, i));
    printf("];\n");
    printf("e2=[");
    for (int i = 0; i < errors2.cols; ++i)
        printf("%f,", errors2.at<float>(0, i));
    printf("];\n");

    double t = findThreshold(errors1, errors2);
    printf("T = %f\n", t);
}

#endif // TEST_HPP