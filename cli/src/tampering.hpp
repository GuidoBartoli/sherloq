#ifndef TAMPERING_HPP
#define TAMPERING_HPP

#include <opencv2/opencv.hpp>
#include "interpolation.h"
#include <opencv2/nonfree/nonfree.hpp>
#include "jpeg.hpp"
#include "utility.hpp"

using namespace std;
using namespace cv;

float angleDistance(float a, float b)
{
    float dist = fabs(a - b);
    if (dist > 180)
        dist = 360 - dist;
    return dist;
}

void cloneDetection(const Mat input, Mat &output)
{
    const int   MAX_NORM    = 135;
    const float MIN_DIST    = 0.01;
    const float MAX_DIST    = 0.05;
    const int   MAX_ROT     = 35;
    const int   MAX_ANGLE   = 4;
    const int   LINE_WIDTH  = 1;

    Mat gray;
    cvtColor(input, gray, CV_BGR2GRAY);
    vector<KeyPoint> keypoints;
    Mat descriptors;
    SIFT sift;
    sift(gray, noArray(), keypoints, descriptors);
    BFMatcher matcher(NORM_L2, true);
    vector< vector<DMatch> > matches;
    matcher.radiusMatch(descriptors, descriptors, matches, MAX_NORM);

    float min_dist = MIN_DIST * min(input.rows, input.cols);
    float max_dist = MAX_DIST * min(input.rows, input.cols);
    // printf("MIN_DIST = %f, MAX_DIST = %f\n", min_dist, max_dist);
    vector< vector<float> > pairs;
    float max_response = numeric_limits<float>::min();
    for (uint i = 0; i < matches.size(); ++i)
    {
        for (uint j = 0; j < matches.at(i).size(); ++j)
        {
            DMatch dmatch = matches.at(i).at(j);
            KeyPoint query = keypoints.at(dmatch.queryIdx);
            KeyPoint train = keypoints.at(dmatch.trainIdx);
            Point p1 = query.pt;
            Point q1 = train.pt;
            double distance = norm(p1 - q1);
            double angle = angleDistance(query.angle, train.angle);
            if (distance > min_dist && angle < MAX_ROT)
            {
                bool found = false;
                for (uint k = 0; k < pairs.size(); ++k)
                {
                    Point p2(pairs.at(k).at(0), pairs.at(k).at(1));
                    Point q2(pairs.at(k).at(2), pairs.at(k).at(3));
                    if ((p1 == p2 && q1 == q2) || (p1 == q2 && q1 == p2))
                    {
                        // printf("FOUND\n");
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    vector<float> kp;
                    kp.push_back(p1.x);
                    kp.push_back(p1.y);
                    kp.push_back(q1.x);
                    kp.push_back(q1.y);
                    kp.push_back(query.size);
                    float response = max(query.response, train.response);
                    kp.push_back(response);
                    if (response > max_response)
                        max_response = response;
                    pairs.push_back(kp);
                    // printf("P = (%d,%d), Q = (%d,%d)\n", p1.x, p1.y, q1.x, q1.y);
                }
            }
        }
    }

    output = input.clone();
    for (uint i = 0; i < pairs.size(); ++i)
    {
        Point p1(pairs.at(i).at(0), pairs.at(i).at(1));
        Point q1(pairs.at(i).at(2), pairs.at(i).at(3));
        double dir1 = atan2(q1.y - p1.y, q1.x - p1.x) * 180 / PI;
        float size = pairs.at(i).at(4);
        float resp = pairs.at(i).at(5) / max_response;
        // printf("P1 = (%d,%d), Q1 = (%d,%d), A1 = %f, ", p1.x, p1.y, q1.x, q1.y, dir1);
        for (uint j = 0; j < pairs.size(); ++j)
        {
            if (j == i)
                continue;
            Point p2(pairs.at(j).at(0), pairs.at(j).at(1));
            Point q2(pairs.at(j).at(2), pairs.at(j).at(3));
            double dir2 = atan2(q2.y - p2.y, q2.x - p2.x) * 180 / PI;
            double p_dist = norm(p1 - p2);
            double q_dist = norm(q1 - q2);
            double dir_diff = angleDistance(dir1, dir2);
            if (dir_diff < MAX_ANGLE && p_dist > min_dist && p_dist < max_dist && q_dist > min_dist && q_dist < max_dist)
            {
                // printf("P2 = (%d,%d), Q2 = (%d,%d), A2 = %f, A12 = %f, D2 = %f", p2.x, p2.y, q2.x, q2.y, dir2, dir_diff, distance2);
                double hue = (dir1 + 180) / 2;
                double sat = resp * 255;
                double val = 255;
                Mat color(1, 1, CV_8UC3);
                color.at<Vec3b>(0, 0) = Vec3b((uchar)hue, (uchar)sat, (uchar)val);
                Mat bgr(1, 1, CV_8UC3);
                cvtColor(color, bgr, CV_HSV2BGR);
                uchar r = bgr.at<Vec3b>(0, 0)[2];
                uchar g = bgr.at<Vec3b>(0, 0)[1];
                uchar b = bgr.at<Vec3b>(0, 0)[0];
                line(output, p1, q1, Scalar(b, g, r), LINE_WIDTH, CV_AA);
                circle(output, p1, size, Scalar(b, g, r), LINE_WIDTH, CV_AA);
                circle(output, q1, size, Scalar(b, g, r), LINE_WIDTH, CV_AA);
                break;
            }
        }
    }
}

void correlation(const Mat input, Mat &plot)
{
    const int   IMG_CHAN  = 1;
    const int   EM_RADIUS = 3;
    const float EM_ERROR  = 0.0001;
    const float EM_STDDEV = 5;
    const int   FLT_SIZE  = 3;
    const int   MAX_ITER  = 20;

    vector<Mat> bgr;
    split(input, bgr);
    Mat y0 = bgr.at(IMG_CHAN);
    y0.convertTo(y0, CV_32F);
    Mat corr(1, y0.cols, CV_32F, Scalar(0));

    for (int r = 0; r < y0.rows; ++r)
    {
        Mat y = y0.row(r).t();
        double min_val, max_val;
        minMaxLoc(y, &min_val, &max_val);
        float p0 = 1 / (max_val - min_val);

        Mat a0(2 * EM_RADIUS + 1, 1, CV_32F, Scalar(1));
        a0.at<float>(EM_RADIUS, 0) = 0;
        a0 /= 2 * EM_RADIUS;
        Mat a1(2 * EM_RADIUS + 1, 1, CV_32F, Scalar(0));

        Mat Y(y.rows, 2 * EM_RADIUS, CV_32F);
        for (int i = 0; i < y.rows; ++i)
        {
            for (int j = 0; j < 2 * EM_RADIUS + 1; ++j)
            {
                if (j == EM_RADIUS)
                    continue;
                Y.at<float>(i, j < EM_RADIUS ? j : j - 1) = y.at<uchar>(i + j, 0);
            }
        }

        int n = 0;
        float s = EM_STDDEV;
        Mat w(y.rows, 1, CV_32F);
        while (norm(a0, a1) > EM_ERROR && n++ < MAX_ITER)
        {
            Mat f;
            filter2D(y, f, CV_32F, a0);
            Mat R = abs(y - f);
            GaussianBlur(R, R, Size(FLT_SIZE, FLT_SIZE), 0);
            pow(R, 2, R);

            Mat W(y.rows, y.rows, CV_32F, Scalar(0));
            float c1 = s * sqrt(2 * PI);
            float c2 = 2 * pow(s, 2);
            float s1 = 0, s2 = 0;
            for (int i = 0; i < R.rows; ++i)
            {
                float Ri = R.at<float>(i, 0);
                float p = (1 / c1) * exp(- Ri / c2);
                float wi = p / (p + p0);
                w.at<float>(i, 0) = wi;
                W.at<float>(i, i) = wi;
                s1 += wi * Ri;
                s2 += wi;
            }

            s = sqrt(s1 / s2);
            a0 = a1.clone();
            Mat a1_0 = (Y.t() * W * Y).inv() * Y.t() * W * y;
            int k = 0;
            for (int i = 0; i < a1_0.cols; ++i)
                if (i != EM_RADIUS)
                    a1.at<float>(i, 0) = a1_0.at<float>(k++, 0);
                else
                    a1.at<float>(i, 0) = 0;
            // printf("  %d) e = %f, s = %f\n", n, norm(a0, a1), s);
        }
        Mat temp;
        transpose(w, temp);
        vconcat(corr, temp, corr);
    }
    corr.rowRange(1, corr.rows).copyTo(corr);
    printMat("CORR", corr, false, true);
    showMat("CORR", corr);
    waitKey(0);
}

void resamplingDetection(const Mat input, Mat &post, Mat &fft, Mat &interp)
{
    const float EM_ERROR  = 0.01;
    const float EM_STDDEV = 5;
    const int   EM_RADIUS = 2;
    const int   FLT_SIZE  = 3;
    const int   MAX_ITER  = 20;
    const int   IMG_SIZE  = 50;
    const int   FFT_LOW   = 140;
    const int   FFT_HIGH  = 150;

    Mat orig;
    cvtColor(input, orig, CV_BGR2GRAY);
    orig.convertTo(orig, CV_32F);
    double min_val, max_val;
    minMaxLoc(orig, &min_val, &max_val);
    float p0 = 1 / (max_val - min_val);
    Mat kernel = (Mat_<float>(FLT_SIZE, FLT_SIZE) << 1, 2, 1, 2, 4, 2, 1, 2, 1);
    kernel /= ((FLT_SIZE + 1) * (FLT_SIZE + 1));
    Mat alpha0 = Mat::ones(2*EM_RADIUS + 1, 2*EM_RADIUS + 1, CV_32F);
    alpha0.at<float>(EM_RADIUS, EM_RADIUS) = 0;
    alpha0 /= (alpha0.rows * alpha0.cols - 1);
    Mat alpha1 = Mat::zeros(2*EM_RADIUS + 1, 2*EM_RADIUS + 1, CV_32F);
    Mat orig_pad(orig, Rect(EM_RADIUS, EM_RADIUS, orig.cols - 2*EM_RADIUS, orig.rows - 2*EM_RADIUS));
    Mat orig0 = orig_pad.clone().reshape(0, orig_pad.rows * orig_pad.cols);

    int block_rows = orig.rows - 2*EM_RADIUS;
    int block_cols = orig.cols - 2*EM_RADIUS;
    Mat resamp(block_rows * block_cols, 1, CV_32F, Scalar(0));
    for (int i = -EM_RADIUS; i <= +EM_RADIUS; ++i)
    {
        for (int j = -EM_RADIUS; j <= +EM_RADIUS; ++j)
        {
            if (i == 0 && j == 0)
                continue;
            Mat block(orig, Rect(j + EM_RADIUS, i + EM_RADIUS, block_cols, block_rows));
            hconcat(resamp, block.clone().reshape(0, block_rows * block_cols), resamp);
        }
    }
    resamp.colRange(1, resamp.cols).copyTo(resamp);
    // resamp.convertTo(resamp, CV_8U);
    // printMat("resamp", resamp);
    // saveMat("resamp.png", resamp);

    int iter = 0;
    float sigma = EM_STDDEV;
    // vector<float> optim;
    // printMat("orig_pad", orig_pad);
    while (norm(alpha0, alpha1) > EM_ERROR && iter++ < MAX_ITER)
    {
        // EXPECTATION
        Mat filt;
        filter2D(orig, filt, CV_32F, alpha1);
        // showMat("FILT", filt);
        Mat resid0 = abs(orig - filt);
        Mat resid(resid0, Rect(EM_RADIUS, EM_RADIUS, resid0.cols - 2*EM_RADIUS, resid0.rows - 2*EM_RADIUS));
        filter2D(resid, resid, CV_32F, kernel);
        // showMat("RESID1", resid);
        // printf("SIZE = %dx%d\n", resid.rows, resid.cols);
        pow(resid, 2, resid);
        // showMat("RESID2", resid);
        Mat cond(resid.rows, resid.cols, CV_32F);
        float c1 = 1 / (sigma * sqrt(2*PI));
        float c2 = 2 * pow(sigma, 2);
        // printf("C1 = %f, C2 = %f\n", c1, c2);
        for (int i = 0; i < cond.rows; ++i)
            for (int j = 0; j < cond.cols; ++j)
                cond.at<float>(i, j) = c1 * (1 / exp(resid.at<float>(i, j) / c2));
        // showMat("COND", cond);
        post = cond / (cond + p0);
        // showMat("POST", post);

        // MAXIMIZATION
        sigma = sqrt(sum(post.mul(resid))[0] / sum(post)[0]) / 2;
        alpha0 = alpha1.clone();
        Mat weights = post.reshape(0, post.rows * post.cols);
        Mat deriv(resamp.cols, resamp.rows, CV_32F);
        for (int i = 0; i < deriv.rows; ++i)
            for (int j = 0; j < deriv.cols; ++j)
                deriv.at<float>(i, j) = resamp.at<float>(j, i) * weights.at<float>(j, 0);
        // showMat("DERIV", deriv, 1/200.0, 10);
        Mat alpha = (deriv * resamp).inv() * deriv * orig0;
        // printMat("ALPHA", alpha.t());

        int k = 0;
        for (int i = 0; i < alpha1.rows; ++i)
            for (int j = 0; j < alpha1.cols; ++j)
                if ((i != alpha1.rows / 2) || (j != alpha1.cols / 2))
                    alpha1.at<float>(i, j) = alpha.at<float>(k++, 0);
                else
                    alpha1.at<float>(i, j) = 0;

        // printMat("ALPHA0", alpha0);
        // showMat("ALPHA0", alpha0, 50, 50);
        // waitKey(50);
        // optim.push_back(norm(alpha0, alpha1));

        // saveMat("cond.png", cond);
        // saveMat("post.png", post);
        // saveMat("resid.png", resid);
        // printMat("POST", post);
        // printMat("COND", cond);
        // printMat("RESID", resid);
        // printMat("ALPHA0", alpha0);
        // printMat("ALPHA", alpha);
        // printMat("ALPHA0", alpha0);

        // printf("resamp(Y) = %dx%d, post(w) = %dx%d, weights(W) = %dx%d\nsignal(y) = %dx%d, alpha0 = %dx%d, orig_pad = %dx%d, alpha = %dx%d, sigma = %f\n", resamp.rows, resamp.cols, post.rows, post.cols, weights.rows, weights.cols, orig.rows, orig.cols, alpha0.rows, alpha0.cols, orig_pad.rows, orig_pad.cols, alpha.rows, alpha.cols, sigma);
        // return;
        // printf("  %d) e = %f, s = %f\n", iter, norm(alpha1 - alpha0), sigma);
    }
    // printMat("ALPHA0", alpha0);
    // printf("optim={");
    // for (uint i = 0; i < optim.size(); ++i)
    //     printf("%f,", optim.at(i));
    // printf("};\n");
    // waitKey();
    normalize(post, post, 0, 1, NORM_MINMAX);

    int M = getOptimalDFTSize( post.rows );
    int N = getOptimalDFTSize( post.cols );
    Mat padded;
    copyMakeBorder(post, padded, 0, M - post.rows, 0, N - post.cols, BORDER_CONSTANT, Scalar::all(0));

    Mat planes[] = {Mat_<float>(padded), Mat::zeros(padded.size(), CV_32F)};
    Mat complexImg;
    merge(planes, 2, complexImg);
    dft(complexImg, complexImg);

    split(complexImg, planes);
    magnitude(planes[0], planes[1], planes[0]);
    Mat mag = planes[0];
    mag += Scalar::all(1);
    log(mag, mag);
    mag = mag(Rect(0, 0, mag.cols & -2, mag.rows & -2));
    int cx = mag.cols/2;
    int cy = mag.rows/2;

    Mat tmp;
    Mat q0(mag, Rect(0, 0, cx, cy));
    Mat q1(mag, Rect(cx, 0, cx, cy));
    Mat q2(mag, Rect(0, cy, cx, cy));
    Mat q3(mag, Rect(cx, cy, cx, cy));
    q0.copyTo(tmp);
    q3.copyTo(q0);
    tmp.copyTo(q3);
    q1.copyTo(tmp);
    q2.copyTo(q1);
    tmp.copyTo(q2);

    normalize(mag, fft, 0, 255, NORM_MINMAX);
    fft.convertTo(fft, CV_8U);
    LUT(fft, buildContrastLUT(FFT_LOW, FFT_HIGH), fft);

    normalize(post, post, 0, 255, NORM_MINMAX);
    post.convertTo(post, CV_8U);

    normalize(alpha1, interp, 0, 255, NORM_MINMAX);
    interp.convertTo(interp, CV_8U);
    resize(interp, interp, Size(), IMG_SIZE, IMG_SIZE, INTER_NEAREST);
}

// ----------------------------------------------------------------------------------- //

Mat histDFT(const Mat hist, bool full = false)
{
    int len = hist.cols;
    Mat complex, planes[] = {Mat_<float>(hist), Mat::zeros(1, len, CV_32F)};
    merge(planes, 2, complex);
    dft(complex, complex);
    split(complex, planes);
    magnitude(planes[0], planes[1], planes[0]);
    Mat dft_mag;
    if (full)
        hconcat(planes[0].colRange(len/2, len), planes[0].colRange(0, len/2), dft_mag);
    else
        dft_mag = planes[0].colRange(0, len/2).clone();
    normalize(dft_mag, dft_mag, 0, 1, NORM_MINMAX);
    return dft_mag;
}

double computeHistErr(const Mat input)
{
    const int    CUT_OFF = 8;
    const double MAX_ERR = 0.185;
    const double SPL_POW = 0.5;

    Mat hist(1, 256, CV_32F, Scalar(0));
    for (int i = 0; i < input.rows; ++i)
        for (int j = 0; j < input.cols; ++j)
            ++hist.at<float>(input.at<uchar>(i, j));

    for (int i = 0; i < hist.cols; ++i)
    {
        if (i <= CUT_OFF)
            hist.at<float>(0, i) *= (1 - cos((PI*(double)i) / CUT_OFF)) / 2;
        else if (i >= 255 - CUT_OFF)
            hist.at<float>(0, i) *= (1 + cos((PI*((double)i+CUT_OFF-255)) / CUT_OFF)) / 2;
    }
    normalize(hist, hist, 0, 1, NORM_MINMAX);

    vector<double> x(256), y(256);
    for (int i = 0; i < 256; ++i)
    {
        x[i] = i;
        y[i] = hist.at<float>(i);
    }

    Mat hist_dft = histDFT(hist, true);
    double en = 0, ed = 0, max_diff = numeric_limits<double>::min();
    for (int i = 0; i < hist.cols; ++i)
    {
        double w = pow(((double)i - 128) / 128, 2);
        en += fabs(w * hist_dft.at<float>(0, i));
        ed += fabs(hist_dft.at<float>(0, i));

        double v = y.at(i);
        x.erase(x.begin() + i);
        y.erase(y.begin() + i);
        alglib::spline1dinterpolant spline;
        alglib::real_1d_array ax, ay;
        ax.setcontent(x.size(), &(x[0]));
        ay.setcontent(y.size(), &(y[0]));
        alglib::spline1dbuildcubic(ax, ay, spline);
        double s = alglib::spline1dcalc(spline, i);
        double diff = fabs(v - s);
        if (diff > max_diff)
            max_diff = diff;
        x.insert(x.begin() + i, i);
        y.insert(y.begin() + i, v);
    }

    double err = en / ed;
    if (err > MAX_ERR)
        err = 1;
    else
        err /= MAX_ERR;

    return (err * pow(max_diff, SPL_POW));
}

double computeChanSim(const Mat input)
{
    const double MAX_SIM = 0.75;

    vector<Mat> bgr, deriv;
    split(input, bgr);
    Mat kx, ky;
    getDerivKernels(kx, ky, 1, 1, 1);
    for (uint c = 0; c < bgr.size(); ++c)
    {
        Mat flt;
        sepFilter2D(bgr.at(c), flt, CV_32F, kx, ky);
        deriv.push_back(flt);
    }

    double s = 0, m = 0;
    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            float b = deriv.at(0).at<float>(i, j);
            float g = deriv.at(1).at<float>(i, j);
            float r = deriv.at(2).at<float>(i, j);
            s += (fabs(g - r) + fabs(g - b) + fabs(r - b));
            m += (fabs(r) + fabs(g) + fabs(b));
        }
    }
    s /= (3 * input.rows * input.cols);
    m /= (3 * input.rows * input.cols);
    double sim = s / m;

    return (sim > MAX_SIM ? 1 : sim / MAX_SIM);
}

void contrastEnhancement(const Mat input, Mat &output)
{
    const int BLOCK_SIZE = min(input.rows, input.cols) / 128 > 10 ? 128 : 64;

    Mat color = padImage(input, BLOCK_SIZE);
    Mat gray;
    cvtColor(color, gray, CV_BGR2GRAY);
    output = Mat(color.rows / BLOCK_SIZE, color.cols / BLOCK_SIZE, CV_32F);
    // Mat output1(input0.rows, input0.cols, CV_32F);
    // Mat output2(input0.rows, input0.cols, CV_32F);
    for (int i = 0; i < color.rows; i += BLOCK_SIZE)
    {
        for (int j = 0; j < color.cols; j += BLOCK_SIZE)
        {
            Mat color_blk(color, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            Mat gray_blk(gray, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            double err = computeHistErr(gray_blk);
            // double spl = computeHistSpl(gray_blk);
            double sim = computeChanSim(color_blk);
            // Mat block1(output1, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            // block1 = err;
            // Mat block2(output2, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            // block2 = sim;
            output.at<float>(i / BLOCK_SIZE, j / BLOCK_SIZE) = err * sim;
        }
    }
    medianBlur(output, output, 3);
    resize(output, output, Size(), BLOCK_SIZE, BLOCK_SIZE, INTER_NEAREST);
    output.convertTo(output, CV_8U, 255);

    // normalize(output1, output1, 0, 255, NORM_MINMAX);
    // normalize(output2, output2, 0, 255, NORM_MINMAX);
    // output1.convertTo(output1, CV_8U);
    // output2.convertTo(output2, CV_8U);
    // imwrite("err.jpg", output1);
    // imwrite("sim.jpg", output2);
}

// ----------------------------------------------------------------------------------- //

void splicingDetection(const Mat input, Mat &output)
{
    Mat gray = padImage(input, DCT_SIZE);
    output = Mat(gray.rows, gray.cols, CV_8UC3);
    for (int i = 0; i < gray.rows; i += DCT_SIZE)
    {
        for (int j = 0; j < gray.cols; j += DCT_SIZE)
        {
            Rect square(j, i, DCT_SIZE, DCT_SIZE);
            int n, l1, l2, q1, q2;
            float p1, p2;
            Mat plot;
            Mat block1(gray, square);
            compressionGhosts(block1, n, l1, p1, l2, p2, q1, q2, plot);
            imwrite("plot.jpg", plot);
            return;
            Mat block2(output, square);
            printf("row = %d, col = %d, n = %d, p1 = %f, p2 = %f\n", i, j, n, p1, p2);
            if (n == 0)
            {
                Mat green(DCT_SIZE, DCT_SIZE, CV_8UC3, Scalar(0, 1, 0));
                green *= (uchar)(p1 / 100 * 255);
                block2 = green;
            }
            else if (n == 1)
            {
                Mat yellow(DCT_SIZE, DCT_SIZE, CV_8UC3, Scalar(0, 1, 1));
                yellow *= (uchar)(p1 / 100 * 255);
                block2 = yellow;
            }
            else
            {
                Mat red(DCT_SIZE, DCT_SIZE, CV_8UC3, Scalar(0, 0, 1));
                red *= (uchar)(p2 / 100 * 255);
                block2 = red;
            }
        }
    }
}

void splicingDetection2(const string filename, Mat &output)
{
    const int NUM_MODES = 15;

    Mat luma_qt, chroma_qt;
    float luma_q, chroma_q;
    int quality;
    if (qualityEstimation(filename, luma_qt, luma_q, chroma_qt, chroma_q, quality) != 0)
        return;
    luma_qt.convertTo(luma_qt, CV_32F);

    Mat gray = padImage(imread(filename, CV_LOAD_IMAGE_GRAYSCALE), DCT_SIZE);
    gray.convertTo(gray, CV_32F);
    Mat coeffs_d(gray.rows, gray.cols, CV_32F, Scalar(0));
    Mat coeffs_q(gray.rows, gray.cols, CV_32F, Scalar(0));
    Mat dct_err(gray.rows, gray.cols, CV_32F, Scalar(0));
    for (int i = 0; i < gray.rows; i += DCT_SIZE)
    {
        for (int j = 0; j < gray.cols; j += DCT_SIZE)
        {
            Rect square(j, i, DCT_SIZE, DCT_SIZE);
            Mat block_g(gray, square);
            Mat coeffs(DCT_SIZE, DCT_SIZE, CV_32F);
            dct(block_g - 128, coeffs);
            Mat block_d(coeffs_d, square);
            block_d = coeffs.clone();
            coeffs /= luma_qt;
            coeffs.convertTo(coeffs, CV_16S);
            coeffs.convertTo(coeffs, CV_32F);
            Mat block_q(coeffs_q, square);
            block_q = coeffs.clone();
            Mat pixels;
            dct(coeffs, pixels, DCT_INVERSE);
            printMat("BLOCK_D", block_d);
            printMat("BLOCK_Q", block_q);
            printMat("PIXELS", pixels);
            return;
        }
    }

    for (int i = 0; i < NUM_MODES; ++i)
    {

    }
}


#endif // TAMPERING_HPP