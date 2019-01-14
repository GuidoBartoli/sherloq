#ifndef LUMINANCE_HPP
#define LUMINANCE_HPP

#include <opencv2/opencv.hpp>

#include "utility.hpp"
//#include "interpolation.h"

using namespace std;
using namespace cv;

void luminanceGradient(const Mat input, Mat &output)
{
    Mat gray;
    cvtColor(input, gray, CV_BGR2GRAY);
    Mat kernel = (Mat_<float>(3, 1) << +1, 0, -1);

    Mat red;
    filter2D(gray, red, CV_32F, kernel);
    normalize(red, red, 0, 255, NORM_MINMAX);
    red.convertTo(red, CV_8U);
    equalizeHist(red, red);

    Mat green;
    filter2D(gray, green, CV_32F, kernel.t());
    normalize(green, green, 0, 255, NORM_MINMAX);
    green.convertTo(green, CV_8U);
    equalizeHist(green, green);

    Mat blue1, blue2;
    Sobel(gray, blue1, CV_32F, 1, 0, CV_SCHARR);
    Sobel(gray, blue2, CV_32F, 0, 1, CV_SCHARR);
    Mat blue = abs(blue1 + blue2);
    normalize(blue, blue, 0, 255, NORM_MINMAX);
    blue.convertTo(blue, CV_8U);
    Mat lut = buildContrastLUT();
    LUT(blue, lut, blue);

    vector<Mat> gradient;
    gradient.push_back(blue);
    gradient.push_back(green);
    gradient.push_back(red);
    merge(gradient, output);
}

void luminanceDetail(const Mat input, Mat &detail)
{
    Mat gray;
    cvtColor(input, gray, CV_BGR2GRAY);
    detail = Mat::zeros(gray.rows, gray.cols, CV_32F);
    for (int i = 1; i < gray.rows - 1; ++i)
    {
        for (int j = 1; j < gray.cols - 1; ++j)
        {
            float l = gray.at<uchar>(i, j - 1);
            float t = gray.at<uchar>(i - 1, j);
            float r = gray.at<uchar>(i, j + 1);
            float b = gray.at<uchar>(i + 1, j);
            float a = (l + t + r + b) / 4;
            float c = gray.at<uchar>(i, j);
            detail.at<float>(i, j) = fabs(c - a);
        }
    }
    normalize(detail, detail, 0, 255, NORM_MINMAX);
    detail.convertTo(detail, CV_8U);
    equalizeHist(detail, detail);
}

void edgeEchoFilter(const Mat input, Mat &edges)
{
    const float SIZE_RATIO = 0.01;
    const int   MAX_SIZE   = 31;
    const int   MIN_SIZE   = 1;

    int size = (int)((input.rows + input.cols) / 2 * SIZE_RATIO);
    size += size % 2 != 0 ? 0 : 1;
    if (size > MAX_SIZE)
        size = MAX_SIZE;
    else if (size < MIN_SIZE)
        size = MIN_SIZE;

    vector<Mat> channels, laplacians;
    split(input, channels);
    Mat lut = buildContrastLUT();
    for (uint i = 0; i < channels.size(); ++i)
    {
        Mat deriv;
        Laplacian(channels.at(i), deriv, CV_32F, size);
        deriv = abs(deriv);
        normalize(deriv, deriv, 0, 255, NORM_MINMAX);
        deriv.convertTo(deriv, CV_8U);
        LUT(deriv, lut, deriv);
        laplacians.push_back(deriv.clone());
    }
    merge(laplacians, edges);
}

// ----------------------------------------------------------------------------------- //

void demosaicAlteration(const Mat input, Mat &output)
{
    output.create(input.rows, input.cols, CV_8UC3);
    output.setTo(0);
    for (int i = 2; i < input.rows - 2; ++i)
    {
        for (int j = 2; j < input.cols - 2; ++j)
        {
            uchar b, g, r;
            b = g = r = 0;
            if ((i % 2 != 0) && (j % 2 != 0))
            {
                // Red pixels are interpolated from distant up/down/left/right
                float r0 = input.at<Vec3b>(i, j)[2];
                float r1 = input.at<Vec3b>(i, j - 2)[2];
                float r2 = input.at<Vec3b>(i - 2, j)[2];
                float r3 = input.at<Vec3b>(i, j + 2)[2];
                float r4 = input.at<Vec3b>(i + 2, j)[2];
                float r_avg = (r1 + r2 + r3 + r4) / 4;
                r = (uchar)(abs(r0 - r_avg));

                // Green pixels are interpolated from close up/down/left/right
                float g0 = input.at<Vec3b>(i, j)[1];
                float g1 = input.at<Vec3b>(i, j - 1)[1];
                float g2 = input.at<Vec3b>(i - 1, j)[1];
                float g3 = input.at<Vec3b>(i, j + 1)[1];
                float g4 = input.at<Vec3b>(i + 1, j)[1];
                float g_avg = (g1 + g2 + g3 + g4) / 4;
                g = (uchar)(abs(g0 - g_avg) * 4);

                // Blue pixels are interpolated from close topleft/topright/bottomright/bottomleft
                float b0 = input.at<Vec3b>(i, j)[0];
                float b1 = input.at<Vec3b>(i - 1, j - 1)[0];
                float b2 = input.at<Vec3b>(i - 1, j + 1)[0];
                float b3 = input.at<Vec3b>(i + 1, j + 1)[0];
                float b4 = input.at<Vec3b>(i + 1, j - 1)[0];
                float b_avg = (b1 + b2 + b3 + b4) / 4;
                b = (uchar)(abs(b0 - b_avg) * 4);
            }
            else if ((i % 2 == 0) && (j % 2 == 0))
            {
                // Blue pixels are interpolated from distant up/down/left/right
                float b0 = input.at<Vec3b>(i, j)[0];
                float b1 = input.at<Vec3b>(i, j - 2)[0];
                float b2 = input.at<Vec3b>(i - 2, j)[0];
                float b3 = input.at<Vec3b>(i, j + 2)[0];
                float b4 = input.at<Vec3b>(i + 2, j)[0];
                float b_avg = (b1 + b2 + b3 + b4) / 4;
                b = (uchar)(abs(b0 - b_avg));

                // Green pixels are interpolated from close up/down/left/right
                float g0 = input.at<Vec3b>(i, j)[1];
                float g1 = input.at<Vec3b>(i, j - 1)[1];
                float g2 = input.at<Vec3b>(i - 1, j)[1];
                float g3 = input.at<Vec3b>(i, j + 1)[1];
                float g4 = input.at<Vec3b>(i + 1, j)[1];
                float g_avg = (g1 + g2 + g3 + g4) / 4;
                g = (uchar)(abs(g0 - g_avg) * 4);

                // Red pixels are interpolated from close topleft/topright/bottomright/bottomleft
                float r0 = input.at<Vec3b>(i, j)[2];
                float r1 = input.at<Vec3b>(i - 1, j - 1)[2];
                float r2 = input.at<Vec3b>(i - 1, j + 1)[2];
                float r3 = input.at<Vec3b>(i + 1, j + 1)[2];
                float r4 = input.at<Vec3b>(i + 1, j - 1)[2];
                float r_avg = (r1 + r2 + r3 + r4) / 4;
                r = (uchar)(abs(r0 - r_avg) * 4);
            }
            else
            {
                // Green pixels are interpolated from distant topleft/topright/bottomright/bottomleft
                float g0 = input.at<Vec3b>(i, j)[1];
                float g1 = input.at<Vec3b>(i - 2, j - 2)[1];
                float g2 = input.at<Vec3b>(i - 2, j + 2)[1];
                float g3 = input.at<Vec3b>(i + 2, j + 2)[1];
                float g4 = input.at<Vec3b>(i + 2, j - 2)[1];
                float g_avg = (g1 + g2 + g3 + g4) / 4;
                g = (uchar)(abs(g0 - g_avg));

                if ((i % 2 != 0) && (j % 2 == 0))
                {
                    // Green pixels among up/down blue pixels and left/right red pixels
                    float b0 = input.at<Vec3b>(i, j)[1];
                    float b1 = input.at<Vec3b>(i - 1, j)[1];
                    float b2 = input.at<Vec3b>(i + 1, j)[1];
                    float b_avg = (b1 + b2) / 2;
                    b = (uchar)(abs(b0 - b_avg) * 2);

                    float r0 = input.at<Vec3b>(i, j)[2];
                    float r1 = input.at<Vec3b>(i, j - 1)[2];
                    float r2 = input.at<Vec3b>(i, j + 1)[2];
                    float r_avg = (r1 + r2) / 2;
                    r = (uchar)(abs(r0 - r_avg) * 2);
                }
                else if ((i % 2 == 0) && (j % 2 != 0))
                {
                    // Green pixels among up/down red pixels and left/right blue pixels
                    float r0 = input.at<Vec3b>(i, j)[2];
                    float r1 = input.at<Vec3b>(i - 1, j)[2];
                    float r2 = input.at<Vec3b>(i + 1, j)[2];
                    float r_avg = (r1 + r2) / 2;
                    r = (uchar)(abs(r0 - r_avg) * 2);

                    float b0 = input.at<Vec3b>(i, j)[0];
                    float b1 = input.at<Vec3b>(i, j - 1)[0];
                    float b2 = input.at<Vec3b>(i, j + 1)[0];
                    float b_avg = (b1 + b2) / 2;
                    b = (uchar)(abs(b0 - b_avg) * 2);
                }
            }
            output.at<Vec3b>(i, j) = Vec3b(b, g, r);
        }
    }
    normalize(output, output, 0, 255, NORM_MINMAX);
}

void dctMap(const Mat input, Mat &output)
{
    const int DCT_SIZE = 8;

    Mat gray;
    cvtColor(padImage(input, DCT_SIZE), gray, CV_BGR2GRAY);
    gray.convertTo(gray, CV_32F);
    output = Mat(gray.rows, gray.cols, CV_32F, Scalar(0));
    for (int i = 0; i < gray.rows; i += DCT_SIZE)
    {
        for (int j = 0; j < gray.cols; j += DCT_SIZE)
        {
            Mat block1(gray, Rect(j, i, DCT_SIZE, DCT_SIZE));
            // block1 -= 128;
            Mat block2(output, Rect(j, i, DCT_SIZE, DCT_SIZE));
            dct(block1, block2);
            block2.at<float>(0, 0) = 0;
            // block2 = abs(block2);
        }
    }
    normalize(output, output, 0, 255, NORM_MINMAX);
    output.convertTo(output, CV_8U);
    equalizeHist(output, output);
}


void cfaInterpolation(const Mat input, Mat &output)
{
    /*
    const int CHANNEL = 1;

    vector<Mat> bgr;
    split(input, bgr);
    Mat channel = bgr.at(CHANNEL);
    channel.convertTo(channel, CV_32F, 1 / 255.0);

    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            if ((i % 2 == 0 && j % 2 == 0) || (i % 2 != 0 && j % 2 != 0))
            {
                vector<double> xv, yv;


                if (i == 0 || j == 0)
                {
                    xv.push_back(j);
                    yv.push_back(i);
                    printf("GRID = 1, ");
                }
                else
                    printf("GRID = 0, ");
                fv.push_back(channel.at<float>(i, j));
                printf("I = %d, J = %d\n", i, j);
            }
        }
    }
    printf("SIZE_X = %d, SIZE_Y = %d, SIZE_Z = %d\n", xv.size(), yv.size(), fv.size());

    // printf("ROWS = %d, COLS = %d\n", input.rows, input.cols);
    // printf("X = %d, Y = %d, F = %d\n", xv.size(), yv.size(), fv.size());

    alglib::real_1d_array x, y, f;
    x.setcontent(xv.size(), &xv.front());
    y.setcontent(yv.size(), &yv.front());
    f.setcontent(fv.size(), &fv.front());
    alglib::spline2dinterpolant spline;
    alglib::spline2dbuildbilinearv(x, xv.size(), y, yv.size(), f, 1, spline);

    output = Mat(input.rows, input.cols, CV_32F, Scalar(0));
    for (int i = 0; i < output.rows; i += 2)
    {
        for (int j = 0; j < output.cols; j += 2)
        {
            double v = alglib::spline2dcalc(spline, j, i);
            output.at<float>(i, j) = fabs(v - channel.at<float>(i, j));
            printf("V = %f, G = %f\n", v, channel.at<float>(i, j));
            // recon.at<float>(i, j) = alglib::spline2dcalc(spline, j, i);
        }
    }
    // output = abs(channel - recon);
    normalize(output, output, 0, 255, NORM_MINMAX);
    output.convertTo(output, CV_8U);
    */
}

#endif // LUMINANCE_HPP
