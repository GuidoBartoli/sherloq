/*
*
* sherlock.cpp
*
* SHERLOCK (a digital image forensic tool)
*
* Author:
* Guido Bartoli <guidobartoli80[at]icloud[dot]com>
*
*/

#include <stdio.h>
#include <sys/stat.h>

#include "general.hpp"
#include "file.hpp"
#include "jpeg.hpp"
#include "colors.hpp"
#include "luminance.hpp"
#include "noise.hpp"
#include "tampering.hpp"
#include "utility.hpp"
#include "test.hpp"


using namespace std;

// ----------------------------------------------------------------------------------- //

// TODO:
//
// Aggiungere una funzione Digest
// Migliorare il filtraggio dei punti individuati da Clone Detection e raggrupparli
// Provare sia BACM che FDF per la doppia compressione locale
// Stima SNR con PCA sui blocchi (varianza = autovalore minore della matrice di covarianza)
// Scomposizione piramidale gaussiana o con wavelet per visualizzare i livelli di dettaglio
// Analisi dell'istogramma per individuare modifiche di contrasto e luminosità
// Estrazioe delle tabelle di Huffman dal JPEG
// Doppia compressione globale con DFT sull'istogramma dei coefficienti DCT
// Analisi della coerenza dell'aberrazione cromatica per rilevare tampering locale
// Analisi Fine-Grained degli artefatti nella CFA
// Aumentare la precisione della conversione in HSV calcolando H manualmente senza modulo
// DCT Map in Amped Authenticate (?)
// Neal Krawetz Color Density Distribution (?)
// Pensare ad una futura GUI: https://ampedsoftware.com/authenticate-samples

// ----------------------------------------------------------------------------------- //

int main(int argc, char *argv[])
{
    const string FILE_EXT  = ".jpg";
    const string FOLDER1   = "RESULTS";
    const string FOLDER_IN = "INPUT";
    const int    GROUP_LEN = 4;
    const int    HORIZ_DIV = 90;

    const string NOT_IMPLEMENTED  = "*not implemented*";
    const string EXIFTOOL_ERROR   = "Error while executing ExifTool!";
    const string THUMB_NOT_FOUND  = "Embedded image not found!";
    const string CLARIFAI_ERROR   = "Unable to contact Clarifai service!";
    const string NOT_A_JPEG_FILE  = "Image is not a JPEG file!";
    const string TABLES_NOT_FOUND = "JPEG quantization tables not found!";
    const string SVM_NOT_FOUND    = "SVM model file not found!";
    const string GEOTAG_NOT_FOUND = "GPS data not found!";

    vector<string> GROUPS;
    vector<string> NAMES;
    vector<string> OPTIONS;

    GROUPS.push_back("GENERAL");
    NAMES.push_back("Image Dimensions");          OPTIONS.push_back("-gd");
    NAMES.push_back("Actual Filesize");           OPTIONS.push_back("-gs");
    NAMES.push_back("Last Modification");         OPTIONS.push_back("-gm");
    NAMES.push_back("Automatic Tagging");         OPTIONS.push_back("-gt");

    GROUPS.push_back("FILE");
    NAMES.push_back("Metadata Extraction");       OPTIONS.push_back("-fm");
    NAMES.push_back("EXIF File Structure");       OPTIONS.push_back("-fs");
    NAMES.push_back("Thumbnail Extraction");      OPTIONS.push_back("-ft");
    NAMES.push_back("Geolocation Data");          OPTIONS.push_back("-fg");

    GROUPS.push_back("JPEG");
    NAMES.push_back("Quality Estimation");        OPTIONS.push_back("-jq");
    NAMES.push_back("Compression Ghosts");        OPTIONS.push_back("-jg");
    NAMES.push_back("Double Compression");        OPTIONS.push_back("-jd");
    NAMES.push_back("Error Level Analysis");      OPTIONS.push_back("-je");

    GROUPS.push_back("COLORS");
    NAMES.push_back("RGB/HSV/Histogram Plots");   OPTIONS.push_back("-ch");
    NAMES.push_back("RGB/YCbCr/Lab/HSV Spaces");  OPTIONS.push_back("-cv");
    NAMES.push_back("RGB PCA Distance/Cross");    OPTIONS.push_back("-cp");
    NAMES.push_back("Min/Max/Avg RGB Values");    OPTIONS.push_back("-cm");

    GROUPS.push_back("LUMINANCE");
    NAMES.push_back("Luminance Gradient");        OPTIONS.push_back("-lg");
    NAMES.push_back("Luminance Detail");          OPTIONS.push_back("-ld");
    NAMES.push_back("Echo Edge Filter");          OPTIONS.push_back("-le");
    NAMES.push_back("DCT Block Map");             OPTIONS.push_back("-lm");

    GROUPS.push_back("NOISE");
    NAMES.push_back("Noise Separation");          OPTIONS.push_back("-ns");
    NAMES.push_back("Min/Max Deviation");         OPTIONS.push_back("-nm");
    NAMES.push_back("LSB Visualization");         OPTIONS.push_back("-nl");
    NAMES.push_back("SNR Consistency");           OPTIONS.push_back("-ns");

    GROUPS.push_back("TAMPERING");
    NAMES.push_back("Contrast Enhancement");      OPTIONS.push_back("-te");
    NAMES.push_back("Clone Region Detection");    OPTIONS.push_back("-tc");
    NAMES.push_back("Resampling Detection");      OPTIONS.push_back("-tr");
    NAMES.push_back("Splicing Detection");        OPTIONS.push_back("-ts");

    // GROUPS.push_back("PROCESS");
    // NAMES.push_back("Levels Adjustment");               OPTIONS.push_back("-pl");
    // NAMES.push_back("HSL Adjustment");                  OPTIONS.push_back("-ph");
    // NAMES.push_back("Histogram Equalization");          OPTIONS.push_back("-pe");
    // NAMES.push_back("Channel Operations");              OPTIONS.push_back("-pc");

    if (argc < 2)
    {
        printf(
            "SHERLOCK v0.27 {author: Guido Bartoli, license: CC BY-SA 4.0}\n"
            "A suite of forensic tools for analyzing and inspecting digital images.\n"
            "\n"
            "Usage: sherlock [ANALYSES] IMAGE\n"
            "\n"
            "Perform ANALYSES (see below) on the provided IMAGE file.\n"
            "If ANALYSES is omitted, all tests are performed.\n"
            "Results are saved in \"%s\" format into \"%s\" folder.\n",
            FILE_EXT.c_str(), FOLDER1.c_str());
        printf("\n");

        int k = 0;
        for (uint i = 0; i < GROUPS.size(); ++i)
        {
            printf("  %s:\n", GROUPS.at(i).c_str());
            for (uint j = 0; j < GROUP_LEN; ++j)
            {
                printf("    %s    %s\n",OPTIONS.at(k).c_str(), NAMES.at(k).c_str());
                ++k;
            }
            printf("    -%c     \"", GROUPS.at(i).c_str()[0]);
            for (uint j = i * GROUP_LEN; j < i * GROUP_LEN + GROUP_LEN; ++j)
            {
                printf("%s", OPTIONS.at(j).c_str());
                if (j < i * GROUP_LEN + GROUP_LEN - 1)
                    printf(" ");
            }
            printf("\"\n");
            if (i != GROUPS.size() - 1)
                printf("\n");
        }
        return 0;
    }

    string filename(argv[argc - 1]);
    if (!isFilenameValid(filename))
    {
        printf("# Filename should not contain spaces!\n");
        return 1;
    }
    if (!isFilePresent(filename))
    {
        printf("# File not found!\n");
        return 2;
    }
    Mat image = imread(filename, CV_LOAD_IMAGE_COLOR);
    if (image.empty())
    {
        printf("# Unable to load image!\n");
        return 3;
    }

    vector<bool> analyses(OPTIONS.size());
    if (argc == 2)
    {
        fill(analyses.begin(), analyses.end(), true);
    }
    else
    {
        fill(analyses.begin(), analyses.end(), false);
        for (int i = 1; i < argc - 1; ++i)
        {
            bool found = false;
            for (uint j = 0; j < GROUPS.size(); ++j)
            {
                char opt[3];
                sprintf(opt, "-%c", GROUPS.at(j).c_str()[0]);
                if (strcmp(argv[i], opt) == 0)
                {
                    fill(analyses.begin() + GROUP_LEN * j, analyses.begin() + GROUP_LEN * (j + 1), true);
                    found = true;
                    break;
                }
            }
            if (found)
                continue;
            uint k = find(OPTIONS.begin(), OPTIONS.end(), string(argv[i])) - OPTIONS.begin();
            if (k != OPTIONS.size())
            {
                analyses[k] = true;
                continue;
            }
            printf("# Unknown option: '%s'\n", argv[i]);
        }
    }
    if (find(analyses.begin(), analyses.end(), true) == analyses.end())
        return 0;

    // ------------------------------------------------------------------------------- //

    deleteFolder(FOLDER1);
    createFolder(FOLDER1);
    for (uint i = 1; i < GROUPS.size(); ++i)
        createFolder(FOLDER1 + "/" + GROUPS.at(i));
    createFolder(FOLDER1 + "/" + FOLDER_IN);
    copyFile(filename, FOLDER1 + "/" + FOLDER_IN);
    string file_path = getFileFromPath(filename);
    string img_format = getImageFormat(filename);

    printf("Analysis started:\n");
    for (int i = 0; i < HORIZ_DIV; ++i)
        printf("=");
    printf("\n");
    printf("> Input filename ........... \"%s\" (%s format)\n", file_path.c_str(), img_format.c_str());
    for (int i = 0; i < HORIZ_DIV; ++i)
        printf("=");
    printf("\n");

    int result;
    string output;
    string folder2;
    clock_t start = clock();
    for (uint i = 0; i < analyses.size(); ++i)
    {
        if (!analyses.at(i))
            continue;
        folder2 = GROUPS.at(i / GROUP_LEN);
        printf("> %s ", NAMES.at(i).c_str());

        if (i == 0)
        {
            int height, width;
            float mpx;
            getImageSize(filename, height, width, mpx);
            printf("......... %dx%d pixels (%.02f Mpx)\n", height, width, mpx);
        }
        else if (i == 1)
        {
            int filesize = getFileSize(filename);
            string human_bytes = bytesToHuman(filesize);
            printf(".......... %d bytes (%s)\n", filesize, human_bytes.c_str());
        }
        else if (i == 2)
        {
            string modification = getFileDate(filename);
            printf("........ %s\n", modification.c_str());
        }
        else if (i == 3)
        {
            vector<string> img_tags = automaticTagging(filename);
            printf("........ ");
            if (!img_tags.empty())
            {
                int k = 0;
                for (uint i = 0; i < img_tags.size(); ++i)
                {
                    string tag = img_tags.at(i);
                    k += tag.size() + 2;
                    if (k > HORIZ_DIV - 29)
                    {
                        printf("\n                             ");
                        k = tag.size() + 2;
                    }
                    printf("%s%s", tag.c_str(), i < img_tags.size() - 1 ? ", " : "\n");
                }
            }
            else
                printf("%s\n", CLARIFAI_ERROR.c_str());
        }
        if (i == 4)
        {
            output = FOLDER1 + "/" + folder2 + "/metadata.txt";
            result = extractMetadata(filename, output);
            printf("...... ");
            if (result != 0)
                printf("%s\n", EXIFTOOL_ERROR.c_str());
            else
                printf("[%s]\n", output.c_str());
        }
        else if (i == 5)
        {
            output = FOLDER1 + "/" + folder2 + "/structure.html";
            result = fileStructure(filename, output);
            printf("...... ");
            if (result != 0)
                printf("%s\n", EXIFTOOL_ERROR.c_str());
            else
                printf("[%s]\n", output.c_str());
        }
        else if (i == 6)
        {
            output = FOLDER1 + "/" + folder2 + "/thumbnail.jpg";
            result = extractThumbnail(filename, output);
            printf("..... ");
            if (result == -1)
                printf("%s\n", EXIFTOOL_ERROR.c_str());
            else if (result == -2)
                printf("%s\n", THUMB_NOT_FOUND.c_str());
            else
                printf("[%s]\n", output.c_str());
        }
        else if (i == 7)
        {
            output = FOLDER1 + "/" + folder2 + "/gpsdata.txt";
            result = gpsData(filename, output);
            printf("......... ");
            if (result == -1)
                printf("%s\n", EXIFTOOL_ERROR.c_str());
            else if (result == -2 || result == -3)
                printf("%s\n", GEOTAG_NOT_FOUND.c_str());
            else
                printf("[%s]\n", output.c_str());
        }
        else if (i == 8)
        {
            Mat luma, chroma;
            float luma_qt, chroma_qt;
            int jpeg_quality;
            result = qualityEstimation(filename, luma, luma_qt, chroma, chroma_qt, jpeg_quality);
            if (result == -1)
                printf(" ...... %s\n", NOT_A_JPEG_FILE.c_str());
            else if (result == -2)
                printf(" ...... %s\n", EXIFTOOL_ERROR.c_str());
            else if (result == -3)
                printf(" ...... %s\n", TABLES_NOT_FOUND.c_str());
            else
            {
                printf("\n");
                printf("  - Luma QT (average level = %.01f%%):", luma_qt);
                printMat("", luma);
                printf("  - Chroma QT (average level = %.01f%%):", chroma_qt);
                printMat("", chroma);
                printf("  - Final JPEG quality = %d%%\n", jpeg_quality);
            }
        }
        else if (i == 9)
        {
            // testLossless();
            printf("....... ");
            int l1, l2, n, q1, q2;
            float p1, p2;
            Mat plot;
            compressionGhosts(image, n, l1, p1, l2, p2, q1, q2, plot);
            if (l1 < 0)
                printf("Uncompressed (P = %.01f%%, Q = [%d%%,%d%%])\n", p1, q1, q2);
            else
            {
                if (l2 < 0)
                    printf("N1: QF = %d%%, P = %.01f%%\n", l1, p1);
                else
                    printf("N%d: (QF1 = %d%%, P1 = %.01f%%) (QF2 = %d%%, P2 = %.01f%%)\n", n, l1, p1, l2, p2);
            }
            output = FOLDER1 + "/" + folder2 + "/ghost_plot" + FILE_EXT;
            imwrite(output, plot);
        }
        else if (i == 10)
        {
            // createTraining();
            // trainSVM();
            // analyzeProbability();
            printf("....... ");
            float probability;
            result = doubleCompression(image, probability);
            if (result == -1)
                printf("%s\n", SVM_NOT_FOUND.c_str());
            else
                printf("Probability = %s%.01f%%\n", probability > 0 ? "+" : "", probability);

            // createTraining(false);
            // trainSVM();
            // int x, y;
            // gridDetection(image, x, y);

            // Mat features;
            // firstDigitFeatures2(image, features);
            // printMat("Features", features);
        }
        else if (i == 11)
        {
            output = FOLDER1 + "/" + folder2 + "/ela" + FILE_EXT;
            Mat ela;
            errorLevelAnalysis(image, ela);
            imwrite(output, ela);
            printf("..... [%s]\n", output.c_str());
        }
        else if (i == 12)
        {
            output = FOLDER1 + "/" + folder2 + "/";
            printf(".. [close plot window(s) to continue]\n");
            Mat lum_hist_img, red_hist_img, green_hist_img, blue_hist_img, comp_hist_img;
            colorPlots(image, lum_hist_img, red_hist_img, green_hist_img, blue_hist_img, comp_hist_img);
            imwrite(output + "hist_lum" + FILE_EXT, lum_hist_img);
            imwrite(output + "hist_red" + FILE_EXT, red_hist_img);
            imwrite(output + "hist_green" + FILE_EXT, green_hist_img);
            imwrite(output + "hist_blue" + FILE_EXT, blue_hist_img);
            imwrite(output + "hist_comp" + FILE_EXT, comp_hist_img);
        }
        else if (i == 13)
        {
            output = FOLDER1 + "/" + folder2 + "/";
            Mat r, g, b, h, s, v, y, cb, cr, l_, a_, b_;
            spaceConversion(image, r, g, b, h, s, v, y, cb, cr, l_, a_, b_);
            imwrite(output + "rgb_r" + FILE_EXT, r);
            imwrite(output + "rgb_g" + FILE_EXT, g);
            imwrite(output + "rgb_b" + FILE_EXT, b);
            imwrite(output + "hsv_h" + FILE_EXT, h);
            imwrite(output + "hsv_s" + FILE_EXT, s);
            imwrite(output + "hsv_v" + FILE_EXT, v);
            imwrite(output + "ycbcr_y" + FILE_EXT, y);
            imwrite(output + "ycbcr_cb" + FILE_EXT, cb);
            imwrite(output + "ycbcr_cr" + FILE_EXT, cr);
            imwrite(output + "lab_l" + FILE_EXT, l_);
            imwrite(output + "lab_a" + FILE_EXT, a_);
            imwrite(output + "lab_b" + FILE_EXT, b_);
            printf(". [%s{rgb,hsv,ycbcr,lab}%s]\n", output.c_str(), FILE_EXT.c_str());
        }
        else if (i == 14)
        {
            output = FOLDER1 + "/" + folder2 + "/pca";
            printf("... [%s{n}%s]\n", output.c_str(), FILE_EXT.c_str());
            vector<Mat> pca;
            colorPCA(image, pca);
            for (uint i = 0; i < pca.size(); ++i)
                imwrite(output + to_string(i + 1) + FILE_EXT, pca.at(i));
        }
        else if (i == 15)
        {
            output = FOLDER1 + "/" + folder2 + "/rgb_";
            Mat rgb_min, rgb_max, rgb_avg;
            minMaxAvgRGB(image, rgb_min, rgb_max, rgb_avg);
            imwrite(output + "min" + FILE_EXT, rgb_min);
            imwrite(output + "max" + FILE_EXT, rgb_max);
            imwrite(output + "avg" + FILE_EXT, rgb_avg);
            printf("... [%s{min,max,avg}%s]\n", output.c_str(), FILE_EXT.c_str());
        }
        else if (i == 16)
        {
            output = FOLDER1 + "/" + folder2 + "/lg" + FILE_EXT;
            Mat lg;
            luminanceGradient(image, lg);
            imwrite(output, lg);
            printf("....... [%s]\n", output.c_str());
        }
        else if (i == 17)
        {
            output = FOLDER1 + "/" + folder2 + "/detail" + FILE_EXT;
            Mat detail;
            luminanceDetail(image, detail);
            imwrite(output, detail);
            printf("......... [%s]\n", output.c_str());
        }
        else if (i == 18)
        {
            output = FOLDER1 + "/" + folder2 + "/edges" + FILE_EXT;
            Mat edges;
            edgeEchoFilter(image, edges);
            imwrite(output, edges);
            printf("......... [%s]\n", output.c_str());
        }
        else if (i == 19)
        {
            // printf("........ %s\n", NOT_IMPLEMENTED.c_str());

            output = FOLDER1 + "/" + folder2 + "/dct_map" + FILE_EXT;
            Mat dct_map;
            dctMap(image, dct_map);
            imwrite(output, dct_map);
            printf("............ [%s]\n", output.c_str());

            // output = FOLDER1 + "/" + folder2 + "/cfa" + FILE_EXT;
            // Mat cfa;
            // cfaInterpolation(image, cfa);
            // imwrite(output, cfa);

            // output = FOLDER1 + "/" + folder2 + "/demosaic" + FILE_EXT;
            // Mat demosaic;
            // demosaicAlteration(image, demosaic);
            // imwrite(output, demosaic);
            // printf("...... [%s]\n", output.c_str());
        }
        else if (i == 20)
        {
            output = FOLDER1 + "/" + folder2 + "/noise" + FILE_EXT;
            Mat noise;
            noiseSeparation(image, noise);
            imwrite(output, noise);
            printf("......... [%s]\n", output.c_str());
        }
        else if (i == 21)
        {
            output = FOLDER1 + "/" + folder2 + "/min_max" + FILE_EXT;
            Mat min_max;
            minMaxDeviation(image, min_max);
            imwrite(output, min_max);
            printf("........ [%s]\n", output.c_str());
        }
        else if (i == 22)
        {
            output = FOLDER1 + "/" + folder2 + "/lsb" + FILE_EXT;
            Mat lsb;
            lsbVisualization(image, lsb);
            imwrite(output, lsb);
            printf("........ [%s]\n", output.c_str());
        }
        else if (i == 23)
        {
            output = FOLDER1 + "/" + folder2 + "/snr" + FILE_EXT;
            Mat snr;
            localSNR(image, snr);
            imwrite(output, snr);
            printf(".......... [%s]\n", output.c_str());
        }
        else if (i == 24)
        {
            output = FOLDER1 + "/" + folder2 + "/contrast" + FILE_EXT;
            Mat contrast;
            contrastEnhancement(image, contrast);
            imwrite(output, contrast);
            printf("..... [%s]\n", output.c_str());
            // testContrast();

            // Mat dct_map;
            // dctMap(image, dct_map);
            // imwrite(FOLDER1 + "/" + folder2 + "/dct_map" + FILE_EXT, dct_map);

            // Mat plot;
            // correlation(image, plot);

            // Mat kurtosis;
            // kurtosisAnalysis(image, kurtosis);
            // imwrite(FOLDER1 + "/" + folder2 + "/kurtosis" + FILE_EXT, kurtosis);

            // int x, y;
            // gridDetection(image, x, y);
        }
        else if (i == 25)
        {
            output = FOLDER1 + "/" + folder2 + "/clone" + FILE_EXT;
            Mat clone;
            cloneDetection(image, clone);
            imwrite(output, clone);
            printf("... [%s]\n", output.c_str());
        }
        else if (i == 26)
        {
            output = FOLDER1 + "/" + folder2 + "/rsmp_";
            Mat post, fft, interp;
            resamplingDetection(image, post, fft, interp);
            imwrite(output + "map" + FILE_EXT, post);
            imwrite(output + "fft" + FILE_EXT, fft);
            imwrite(output + "int" + FILE_EXT, interp);
            printf("..... [%s{map,fft,int}%s]\n", output.c_str(), FILE_EXT.c_str());
        }
        else if (i == 27)
        {
            printf("....... %s\n", NOT_IMPLEMENTED.c_str());

            // output = FOLDER1 + "/" + folder2 + "/splicing" + FILE_EXT;
            // Mat splicing;
            // splicingDetection(image, splicing);
            // imwrite(output, splicing);
            // printf("....... [%s]\n", output.c_str());
        }
    }
    for (int i = 0; i < HORIZ_DIV; ++i)
        printf("=");
    printf("\n");
    printf("Processing done! (");
    double t = (double)(clock() - start) / CLOCKS_PER_SEC;
    if (t < 1)
        printf("%d ms", (int)(t * 1000));
    else
        printf("%.03f s", t);
    printf(")\n");

    return 0;
}
