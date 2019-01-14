#ifndef COLORS_HPP
#define COLORS_HPP

#include <opencv2/opencv.hpp>
#include <opencv2/viz/viz3d.hpp>

using namespace std;
using namespace cv;

int spaceConversion(const Mat input,
                    Mat &r, Mat &g, Mat &b,
                    Mat &h, Mat &s, Mat &v,
                    Mat &y, Mat &cb, Mat &cr,
                    Mat &l_, Mat &a_, Mat &b_)
{
    vector<Mat> channels;
    split(input, channels);
    r = channels.at(2).clone();
    g = channels.at(1).clone();
    b = channels.at(0).clone();

    Mat input0;
    input.convertTo(input0, CV_32FC3);
    input0 /= 255.0;

    Mat hsv;
    cvtColor(input0, hsv, CV_BGR2HSV);
    split(hsv, channels);
    h = channels.at(0).clone();
    for (int i = 0; i < h.rows; ++i)
    {
        for (int j = 0; j < h.cols; ++j)
        {
            float val = h.at<float>(i, j);
            h.at<float>(i, j) = val > 180 ? 360 - val : val;
        }
    }
    normalize(h, h, 0, 255, NORM_MINMAX);
    h.convertTo(h, CV_8U);
    s = channels.at(1).clone();
    s.convertTo(s, CV_8U, 255);
    v = channels.at(2).clone();
    v.convertTo(v, CV_8U, 255);

    Mat ycbcr;
    cvtColor(input0, ycbcr, CV_BGR2YCrCb);
    split(ycbcr, channels);
    y = channels.at(0).clone();
    y.convertTo(y, CV_8U, 255);
    cr = channels.at(1).clone();
    normalize(cr, cr, 0, 255, NORM_MINMAX);
    cr.convertTo(cr, CV_8U);
    cb = channels.at(2).clone();
    normalize(cb, cb, 0, 255, NORM_MINMAX);
    cb.convertTo(cb, CV_8U);

    Mat lab;
    cvtColor(input0, lab, CV_BGR2Lab);
    split(lab, channels);
    l_ = channels.at(0).clone();
    l_.convertTo(l_, CV_8U, 255 / 100.0);
    a_ = channels.at(1).clone();
    normalize(a_, a_, 0, 255, NORM_MINMAX);
    a_.convertTo(a_, CV_8U);
    b_ = channels.at(2).clone();
    normalize(b_, b_, 0, 255, NORM_MINMAX);
    b_.convertTo(b_, CV_8U);

    return 0;

    /*
    r = Mat(input.rows, input.cols, CV_8U);
    g = Mat(input.rows, input.cols, CV_8U);
    b = Mat(input.rows, input.cols, CV_8U);
    h = Mat(input.rows, input.cols, CV_32F);
    s = Mat(input.rows, input.cols, CV_32F);
    v = Mat(input.rows, input.cols, CV_32F);
    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            Vec3b bgr = input.at<Vec3b>(i, j);
            r.at<uchar>(i, j) = bgr[2];
            g.at<uchar>(i, j) = bgr[1];
            b.at<uchar>(i, j) = bgr[0];
            float r0 = bgr[2] / 255.0;
            float g0 = bgr[1] / 255.0;
            float b0 = bgr[0] / 255.0;
            float v0 = max(max(r0, g0), b0);
            float m = min(min(r0, g0), b0);
            v.at<float>(i, j) = v0;
            s.at<float>(i, j) = v0 != 0 ? (v0 - m) / v0 : 0;
            if (v0 == r0)
                h.at<float>(i, j) = (60 * (g0 - b0)) / (v0 - m);
            else if (v0 == g0)
                h.at<float>(i, j) = 120 + (60 * (b0 - r0)) / (v0 - m);
            else if (v0 == b0)
                h.at<float>(i, j) = 240 + (60 * (r0 - g0)) / (v0 - m);
        }
    }

    Mat hsv;
    cvtColor(input0, hsv, CV_BGR2HSV);
    split(hsv, channels);
    h = channels.at(0).clone();
    for (int i = 0; i < h.rows; ++i)
    {
        for (int j = 0; j < h.cols; ++j)
        {
            float val = h.at<float>(i, j);
            h.at<float>(i, j) = ((val > 180 ? 360 - val : val) / 180);
        }
    }
    s = channels.at(1).clone();
    v = channels.at(2).clone();

    normalize(h, h, 0, 255, NORM_MINMAX);
    h.convertTo(h, CV_8U);
    normalize(s, s, 0, 255, NORM_MINMAX);
    s.convertTo(s, CV_8U);
    normalize(v, v, 0, 255, NORM_MINMAX);
    v.convertTo(v, CV_8U);

    Mat ycbcr;
    cvtColor(input0, ycbcr, CV_BGR2YCrCb);
    split(ycbcr, channels);
    y = channels.at(0).clone();
    cb = channels.at(1).clone();
    normalize(cb, cb, 0, 255, NORM_MINMAX);
    cr = channels.at(2).clone();
    normalize(cr, cr, 0, 255, NORM_MINMAX);

    Mat lab;
    cvtColor(input0, lab, CV_BGR2Lab);
    split(lab, channels);
    l_ = channels.at(0).clone();
    a_ = channels.at(1).clone();
    normalize(a_, a_, 0, 255, NORM_MINMAX);
    b_ = channels.at(2).clone();
    normalize(b_, b_, 0, 255, NORM_MINMAX);

    return 0;
    */
}

void minMaxAvgRGB(const Mat input, Mat &rgb_min, Mat &rgb_max, Mat &rgb_avg)
{
    rgb_min.create(input.rows, input.cols, CV_8UC3);
    rgb_min.setTo(0);
    rgb_max.create(input.rows, input.cols, CV_8UC3);
    rgb_max.setTo(0);
    rgb_avg.create(input.rows, input.cols, CV_8UC3);
    rgb_avg.setTo(0);

    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            Vec3b bgr = input.at<Vec3b>(i, j);
            uchar b = bgr[0];
            uchar g = bgr[1];
            uchar r = bgr[2];

            if (b > g && b > r)
                rgb_max.at<Vec3b>(i, j) = Vec3b(b, 0, 0);
            else if (g > r && g > b)
                rgb_max.at<Vec3b>(i, j) = Vec3b(0, g, 0);
            else if (r > b && r > g)
                rgb_max.at<Vec3b>(i, j) = Vec3b(0, 0, r);

            if (b < g && b < r)
                rgb_min.at<Vec3b>(i, j) = Vec3b(b, 0, 0);
            else if (g < r && g < b)
                rgb_min.at<Vec3b>(i, j) = Vec3b(0, g, 0);
            else if (r < b && r < g)
                rgb_min.at<Vec3b>(i, j) = Vec3b(0, 0, r);

            if ((b > r && b < g) || (b > g && b < r))
                rgb_avg.at<Vec3b>(i, j) = Vec3b(b, 0, 0);
            else if ((g > r && g < b) || (g > b && g < r))
                rgb_avg.at<Vec3b>(i, j) = Vec3b(0, g, 0);
            else if ((r > b && r < g) || (r > g && r < b))
                rgb_avg.at<Vec3b>(i, j) = Vec3b(0, 0, r);
        }
    }
}

void colorPCA(const Mat input, vector<Mat> &output)
{
    const int DIMS = 3;

    Mat rgb_data(input.rows * input.cols, DIMS, CV_32F);
    MatConstIterator_<Vec3b> input_it = input.begin<Vec3b>();
    MatConstIterator_<Vec3b> input_end = input.end<Vec3b>();
    int i = 0;
    while (input_it != input_end)
    {
        for (int j = 0; j < DIMS; ++j)
            rgb_data.at<float>(i, j) = (*input_it)[j];
        ++input_it;
        ++i;
    }
    PCA pca(rgb_data, noArray(), CV_PCA_DATA_AS_ROW, DIMS);

    output.clear();
    Point3f q(pca.mean.at<float>(0), pca.mean.at<float>(1), pca.mean.at<float>(2));
    for (int d = 0; d < DIMS; ++d)
    {
        Point3f v = Point3f(pca.eigenvectors.at<float>(d, 0),
                            pca.eigenvectors.at<float>(d, 1),
                            pca.eigenvectors.at<float>(d, 2));
        Point3f r = q + v * pca.eigenvalues.at<float>(d);
        Mat distance(input.rows, input.cols, CV_32F);
        Mat cross(input.rows, input.cols, CV_32FC3);
        MatIterator_<float> distance_it = distance.begin<float>();
        MatIterator_<Vec3f> cross_it = cross.begin<Vec3f>();
        input_it = input.begin<Vec3b>();
        while (input_it != input_end)
        {
            Point3f p((*input_it)[0], (*input_it)[1], (*input_it)[2]);
            *distance_it = norm((p - q).cross(p - r)) / norm(r - q);
            *cross_it = p.cross(v);
            ++input_it;
            ++distance_it;
            ++cross_it;
        }

        vector<Mat> channels;
        split(cross, channels);
        for (int c = 0; c < 3; ++c)
        {
            Mat channel;
            normalize(channels.at(c), channel, 0, 255, NORM_MINMAX);
            channel.convertTo(channel, CV_8U);
            channels[c] = channel;
        }

        normalize(distance, distance, 0, 255, NORM_MINMAX);
        distance.convertTo(distance, CV_8U);
        output.push_back(distance);

        Mat distance_eq;
        equalizeHist(distance, distance_eq);
        output.push_back(distance_eq);

        Mat cross_img;
        merge(channels, cross_img);
        output.push_back(cross_img);
    }
}

// ----------------------------------------------------------------------------------- //

void colorPlots(const Mat input, Mat &lum_hist_img, Mat &red_hist_img, Mat &green_hist_img, Mat &blue_hist_img, Mat &comp_hist_img)
{
    const int   HIST_HEIGHT = 128;
    const int   HIST_GRAD   = 8;
    const int   HIST_RESIZE = 2;

    const int   GRID_SIZE   = 20;
    const float GRID_SPACE  = 0.1;
    const int   LEGEND_PAD  = 5;
    const int   LEGEND_SIZE = 14;
    const float LABEL_SIZE  = 0.075;
    const float LABEL_DIST  = 1.05;
    const int   POINT_SIZE  = 3;

    Mat lum_hist(1, 256, CV_32S, Scalar(0));
    Mat red_hist = lum_hist.clone();
    Mat green_hist = lum_hist.clone();
    Mat blue_hist = lum_hist.clone();
    lum_hist_img = Mat(HIST_HEIGHT, 256, CV_8UC3, Scalar(0, 0, 0));
    red_hist_img = lum_hist_img.clone();
    green_hist_img = lum_hist_img.clone();
    blue_hist_img = lum_hist_img.clone();

    Mat rgb_points(input.rows * input.cols, 1, CV_32FC3);
    Mat hsv_points(input.rows * input.cols, 1, CV_32FC3);
    Mat rgb_colors(input.rows * input.cols, 1, CV_8UC3);
    Mat hsv, input0;
    input.convertTo(input0, CV_32FC3);
    cvtColor(input0 / 255, hsv, CV_BGR2HSV);

    int k = 0;
    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            uchar b = input.at<Vec3b>(i, j)[0];
            uchar g = input.at<Vec3b>(i, j)[1];
            uchar r = input.at<Vec3b>(i, j)[2];
            float h = hsv.at<Vec3f>(i, j)[0];
            float s = hsv.at<Vec3f>(i, j)[1];
            float v = hsv.at<Vec3f>(i, j)[2];

            hsv_points.at<Vec3f>(k, 0) = Vec3f(h / 360, s, v);
            rgb_points.at<Vec3f>(k, 0) = Vec3f(r / 255.0, g / 255.0, b / 255.0);
            rgb_colors.at<Vec3b>(k, 0) = Vec3b(b, g, r);
            ++k;

            ++lum_hist.at<int>(0, (int)round(v * 255));
            ++red_hist.at<int>(0, r);
            ++green_hist.at<int>(0, g);
            ++blue_hist.at<int>(0, b);
        }
    }

    lum_hist.convertTo(lum_hist, CV_32F);
    normalize(lum_hist, lum_hist, 0, 1, NORM_MINMAX);
    red_hist.convertTo(red_hist, CV_32F);
    normalize(red_hist, red_hist, 0, 1, NORM_MINMAX);
    green_hist.convertTo(green_hist, CV_32F);
    normalize(green_hist, green_hist, 0, 1, NORM_MINMAX);
    blue_hist.convertTo(blue_hist, CV_32F);
    normalize(blue_hist, blue_hist, 0, 1, NORM_MINMAX);

    Mat lum_grad(1, 256, CV_8UC3);
    Mat red_grad = lum_grad.clone();
    Mat green_grad = lum_grad.clone();
    Mat blue_grad = lum_grad.clone();
    for (int i = 0; i < 256; ++i)
    {
        int height = (int)round(lum_hist.at<float>(0, i) * HIST_HEIGHT);
        if (height > 0)
            line(lum_hist_img, Point(i, HIST_HEIGHT - 1), Point(i, (HIST_HEIGHT - 1) - height), Scalar(255, 255, 255));
        height = (int)round(red_hist.at<float>(0, i) * HIST_HEIGHT);
        if (height > 0)
            line(red_hist_img, Point(i, HIST_HEIGHT - 1), Point(i, (HIST_HEIGHT - 1) - height), Scalar(0, 0, 255));
        height = (int)round(green_hist.at<float>(0, i) * HIST_HEIGHT);
        if (height > 0)
            line(green_hist_img, Point(i, HIST_HEIGHT - 1), Point(i, (HIST_HEIGHT - 1) - height), Scalar(0, 255, 0));
        height = (int)round(blue_hist.at<float>(0, i) * HIST_HEIGHT);
        if (height > 0)
            line(blue_hist_img, Point(i, HIST_HEIGHT - 1), Point(i, (HIST_HEIGHT - 1) - height), Scalar(255, 0, 0));
        lum_grad.at<Vec3b>(0, i) = Vec3b(i, i, i);
        red_grad.at<Vec3b>(0, i) = Vec3b(0, 0, i);
        green_grad.at<Vec3b>(0, i) = Vec3b(0, i, 0);
        blue_grad.at<Vec3b>(0, i) = Vec3b(i, 0, 0);
    }
    resize(lum_grad, lum_grad, Size(), 1, HIST_GRAD, INTER_NEAREST);
    resize(red_grad, red_grad, Size(), 1, HIST_GRAD, INTER_NEAREST);
    resize(green_grad, green_grad, Size(), 1, HIST_GRAD, INTER_NEAREST);
    resize(blue_grad, blue_grad, Size(), 1, HIST_GRAD, INTER_NEAREST);
    Mat zero = Mat::zeros(1, 256, CV_8UC3);
    vconcat(lum_hist_img, zero, lum_hist_img);
    vconcat(lum_hist_img, lum_grad, lum_hist_img);
    vconcat(red_hist_img, zero, red_hist_img);
    vconcat(red_hist_img, red_grad, red_hist_img);
    vconcat(green_hist_img, zero, green_hist_img);
    vconcat(green_hist_img, green_grad, green_hist_img);
    vconcat(blue_hist_img, zero, blue_hist_img);
    vconcat(blue_hist_img, blue_grad, blue_hist_img);

    resize(lum_hist_img, lum_hist_img, Size(), HIST_RESIZE, HIST_RESIZE, INTER_NEAREST);
    resize(red_hist_img, red_hist_img, Size(), HIST_RESIZE, HIST_RESIZE, INTER_NEAREST);
    resize(green_hist_img, green_hist_img, Size(), HIST_RESIZE, HIST_RESIZE, INTER_NEAREST);
    resize(blue_hist_img, blue_hist_img, Size(), HIST_RESIZE, HIST_RESIZE, INTER_NEAREST);
    comp_hist_img = red_hist_img + green_hist_img + blue_hist_img;

    // --------------------------------------------------------------------- //

    viz::WGrid grid(Vec2i(GRID_SIZE, GRID_SIZE), Vec2d::all(GRID_SPACE), viz::Color::gray());
    viz::WText legend(
        "Mouse:\n"
        "[L] rotate\n"
        "[R] zoom\n"
        "[M] pan\n"
        "\n"
        "Keyboard:\n"
        "+/- point size\n"
        "[J] save snapshot\n"
        "[R] reset viewpoint\n"
        "[O] change projection\n"
        "[K] export to .obj\n"
        "[E] close window\n",
        Point2i(LEGEND_PAD, LEGEND_PAD), LEGEND_SIZE);

    viz::Viz3d rgb_plot("RGB");
    rgb_plot.showWidget("rgb_origin", viz::WCoordinateSystem());
    rgb_plot.showWidget("rgb_points", viz::WCloud(rgb_points, rgb_colors));
    rgb_plot.setRenderingProperty("rgb_points", viz::POINT_SIZE, POINT_SIZE);
    rgb_plot.showWidget("rgb_grid", grid);
    rgb_plot.showWidget("rgb_legend", legend);
    rgb_plot.spinOnce();
    rgb_plot.showWidget("red_label", viz::WText3D("R", Point3f(LABEL_DIST, 0, 0), LABEL_SIZE));
    rgb_plot.showWidget("green_label", viz::WText3D("G", Point3f(0, LABEL_DIST, 0), LABEL_SIZE));
    rgb_plot.showWidget("blue_label", viz::WText3D("B", Point3f(0, 0, LABEL_DIST), LABEL_SIZE));
    rgb_plot.spinOnce();

    viz::Viz3d hsv_plot("HSV");
    hsv_plot.showWidget("hsv_origin", viz::WCoordinateSystem());
    hsv_plot.showWidget("hsv_points", viz::WCloud(hsv_points, rgb_colors));
    hsv_plot.setRenderingProperty("hsv_points", viz::POINT_SIZE, POINT_SIZE);
    hsv_plot.showWidget("hsv_grid", grid);
    hsv_plot.showWidget("hsv_legend", legend);
    hsv_plot.spinOnce();
    hsv_plot.showWidget("hue_label", viz::WText3D("H", Point3f(LABEL_DIST, 0, 0), LABEL_SIZE));
    hsv_plot.showWidget("sat_label", viz::WText3D("S", Point3f(0, LABEL_DIST, 0), LABEL_SIZE));
    hsv_plot.showWidget("val_label", viz::WText3D("V", Point3f(0, 0, LABEL_DIST), LABEL_SIZE));
    hsv_plot.spin();
}

// ----------------------------------------------------------------------------------- //

void colorDensityDistribution(const Mat input, Mat &output)
{
    const int DIMS = 3;
    const int COLORS = 3;
    const int ITERATIONS = 50;
    const int EPSILON = 0.1;

    Mat rgb_data(input.rows * input.cols, DIMS, CV_32F);
    int r = 0;
    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            for (int d = 0; d < DIMS; ++d)
                rgb_data.at<float>(r, d) = input.at<Vec3b>(i, j)[d];
            ++r;
        }
    }
    Mat labels, centroids;
    TermCriteria criteria(TermCriteria::COUNT + TermCriteria::EPS, ITERATIONS, EPSILON);
    kmeans(rgb_data, COLORS, labels, criteria, 1, KMEANS_PP_CENTERS, centroids);

    output.create(input.rows, input.cols, CV_32FC3);
    int k = 0;
    for (int i = 0; i < output.rows; ++i)
    {
        for (int j = 0; j < output.cols; ++j)
        {
            for (int h = 0; h < COLORS; ++h)
            {
                Vec3f centroid(centroids.at<float>(h, 0),
                               centroids.at<float>(h, 1),
                               centroids.at<float>(h, 2));
                output.at<Vec3f>(i, j)[h] = norm(input.at<Vec3f>(i, j), centroid);
            }
            // int label = labels.at<int>(k, 0);
            // Vec3f centroids
            // uchar b = (uchar)centroids.at<float>(label, 0);
            // uchar g = (uchar)centroids.at<float>(label, 1);
            // uchar r = (uchar)centroids.at<float>(label, 2);
            // output.at<Vec3b>(i, j) = Vec3b(b, g, r);
            ++k;
        }
    }
    // output /= (255 * sqrt(3.0));
    // output *= 255;
    output.convertTo(output, CV_8UC3);
}

void colorDensityDistribution2(const Mat input, Mat &output)
{
    output.create(input.rows, input.cols, CV_32FC3);
    for (int i = 0; i < input.rows; ++i)
    {
        for (int j = 0; j < input.cols; ++j)
        {
            Vec3f bgr(input.at<Vec3b>(i, j)[0],
                      input.at<Vec3b>(i, j)[1],
                      input.at<Vec3b>(i, j)[2]);
            float b = pow(norm(bgr, Vec3f(255, 0, 0)), 2);
            float g = pow(norm(bgr, Vec3f(0, 255, 0)), 2);
            float r = pow(norm(bgr, Vec3f(0, 0, 255)), 2);
            output.at<Vec3f>(i, j) = Vec3f(b, g, r);
        }
    }
    normalize(output, output, 0, 255, NORM_MINMAX);
    output.convertTo(output, CV_8UC3);
}

void highContrastRGB(const Mat input, Mat &output)
{
    const uchar LOW_VALUE = 100;
    const uchar HIGH_VALUE = 145;

    Mat lut = buildContrastLUT(LOW_VALUE, HIGH_VALUE);
    LUT(input, lut, output);
}

#endif // COLORS_HPP