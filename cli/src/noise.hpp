#ifndef NOISE_HPP
#define NOISE_HPP

#include <opencv2/opencv.hpp>
#include "utility.hpp"

using namespace std;
using namespace cv;

void noiseSeparation(const Mat input, Mat &output)
{
    const int   BLOCK_SIZE = 3;
    const int   BLOCK_PAD  = BLOCK_SIZE / 2;
    const float MAX_STDDEV = 135;

    Mat gray, denoised, noise;
    cvtColor(input, gray, CV_BGR2GRAY);
    medianBlur(gray, denoised, BLOCK_SIZE);
    gray.convertTo(gray, CV_32F);
    denoised.convertTo(denoised, CV_32F);
    noise = abs(gray - denoised);

    for (int i = BLOCK_PAD; i < gray.rows - BLOCK_PAD; ++i)
    {
        for (int j = BLOCK_PAD; j < gray.cols - BLOCK_PAD; ++j)
        {
            Mat block(gray, Rect(j - BLOCK_PAD, i - BLOCK_PAD, BLOCK_SIZE, BLOCK_SIZE));
            Scalar mean_sclr, stddev_sclr;
            meanStdDev(block, mean_sclr, stddev_sclr);
            float stddev = stddev_sclr[0] / MAX_STDDEV;
            if (stddev != 0)
                noise.at<float>(i, j) /= stddev;
            else
                noise.at<float>(i, j) = 0;
        }
    }

    GaussianBlur(noise, noise, Size(BLOCK_SIZE, BLOCK_SIZE), 0);
    normalize(noise, noise, 0, 255, NORM_MINMAX);
    noise.convertTo(output, CV_8U);
}

void minMaxDeviation(const Mat input, Mat &output)
{
    const Vec3b MAX_COLOR(0, 0, 255);
    const Vec3b MIN_COLOR(0, 255, 0);

    output = Mat::zeros(input.rows, input.cols, CV_8UC3);
    for (int i = 1; i < input.rows - 1; ++i)
    {
        for (int j = 1; j < input.cols - 1; ++j)
        {
            float norm0 = norm(input.at<Vec3b>(i, j));
            float norm1 = norm(input.at<Vec3b>(i, j - 1));
            float norm2 = norm(input.at<Vec3b>(i - 1, j));
            float norm3 = norm(input.at<Vec3b>(i, j + 1));
            float norm4 = norm(input.at<Vec3b>(i + 1, j));
            float max_norm = max(max(max(norm1, norm2), norm3), norm4);
            float min_norm = min(min(min(norm1, norm2), norm3), norm4);
            if (norm0 > max_norm)
                output.at<Vec3b>(i, j) = MAX_COLOR;
            else if (norm0 < min_norm)
                output.at<Vec3b>(i, j) = MIN_COLOR;
        }
    }
}

void localSNR(const Mat input, Mat &output)
{
    const int BLOCK_SIZE = 8;

    Mat min_max;
    minMaxDeviation(input, min_max);
    output.create(input.rows, input.cols, CV_32F);
    output.setTo(255);
    Mat mask(input.rows, input.cols, CV_8U, Scalar(0));
    for (int i = 0; i < min_max.rows - BLOCK_SIZE; i += BLOCK_SIZE)
    {
        for (int j = 0; j < min_max.cols - BLOCK_SIZE; j += BLOCK_SIZE)
        {
            Mat block1(min_max, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            Scalar avg, stddev;
            meanStdDev(block1, avg, stddev);
            Mat block2(output, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            block2 = stddev[1] + stddev[2];
            Mat block3(mask, Rect(j, i, BLOCK_SIZE, BLOCK_SIZE));
            block3 = 1;
        }
    }
    normalize(output, output, 0, 255, NORM_MINMAX, -1, mask);
    output = 255 - output;
    output.convertTo(output, CV_8U);
}

void lsbVisualization(const Mat input, Mat &output)
{
    output.create(input.rows, input.cols, CV_8UC3);
    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            Vec3b bgr = input.at<Vec3b>(i, j);
            int b = bgr[0];
            int g = bgr[1];
            int r = bgr[2];
            output.at<Vec3b>(i, j) =
                    Vec3b(b % 2 == 0 ? 0 : 255, g % 2 == 0 ? 0 : 255, r % 2 == 0 ? 0 : 255);
        }
    }
}

#endif // NOISE_HPP