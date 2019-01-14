#ifndef FILE_HPP
#define FILE_HPP

#include <opencv2/opencv.hpp>
#include "utility.hpp"

using namespace std;
using namespace cv;

int extractMetadata(const string image_file, const string data_file)
{
    return runExifTool("-g --GPS* --preview:all", image_file, data_file);
}

int fileStructure(const string image_file, const string struct_file)
{
    return runExifTool("-htmlDump0", image_file, struct_file);
}

int extractThumbnail(const string image_file, const string thumb_file)
{
    // exiftool -a -b -W FOLDERNAME/%f_%t%-c.%s -preview:all FILE
    int code = runExifTool("-b -PreviewImage", image_file, thumb_file);
    if (code != 0)
        return -1;
    if (checkImage(thumb_file))
        return 0;
    runExifTool("-b -ThumbnailImage", image_file, thumb_file);
    if (checkImage(thumb_file))
        return 0;
    runExifTool("-b -preview:all", image_file, thumb_file);
    if (checkImage(thumb_file))
        return 0;
    deleteFile(thumb_file);
    return -2;
}

int gpsData(const string image_file, const string output_file)
{
    const string TEMP_FILE     = "temp.txt";
    const int    BUFF_LENGTH   = 256;
    const char   DELIMITERS[6] = {' ', ':', ',', '\'', '\"', '\n'};

    int code = runExifTool("-GPSPosition", image_file, TEMP_FILE);
    if (code != 0)
    {
        deleteFile(TEMP_FILE);
        return -1;
    }
    FILE* file = fopen(TEMP_FILE.c_str(), "r");
    char line[BUFF_LENGTH];
    char* ret = fgets(line, sizeof(line), file);
    fclose(file);
    deleteFile(TEMP_FILE);
    if (ret == NULL)
        return -2;
    // printf("\n%s\n", line);

    double lat = 0, lon = 0;
    char* token = strtok(line, DELIMITERS);
    int i = 0;
    char lat_dir = 0;
    char lon_dir = 0;
    while (token != NULL)
    {
        string str(token);
        float val = atof(token);
        // printf("%d) <%s>\n", i, str.c_str());
        if (i == 2)
            lat += val;
        else if (i == 4)
            lat += val / 60;
        else if (i == 5)
            lat += val / 3600;
        else if (i == 6)
        {
            if (str == "S")
            {
                lat *= -1;
                lat_dir = 'S';
            }
            else if (str == "N")
                lat_dir = 'N';
            else
                return -3;
        }
        else if (i == 7)
            lon += val;
        else if (i == 9)
            lon += val / 60;
        else if (i == 10)
            lon += val / 3600;
        else if (i == 11)
        {
            if (str == "W")
            {
                lon *= -1;
                lon_dir = 'W';
            }
            else if (str == "E")
                lon_dir = 'E';
            else
                return -3;
        }
        token = strtok(NULL, DELIMITERS);
        ++i;
    }
    if (lat == 0 || lon == 0)
        return -3;

    double lat_d = (int)fabs(lat);
    double lat_m = (int)((fabs(lat) - lat_d) * 60);
    double lat_s = (fabs(lat) - lat_d - lat_m / 60) * 3600;
    double lon_d = (int)fabs(lon);
    double lon_m = (int)((fabs(lon) - lon_d) * 60);
    double lon_s = (fabs(lon) - lon_d - lon_m / 60) * 3600;
    string maps_url = "Google Map URL: https://maps.google.com/maps?ll=" + to_string(lat) + "," + to_string(lon) + "&t=h&q=" + to_string((int)lat_d) + "%C2%B0" + to_string((int)lat_m) + "%27" + to_string(lat_s) + "%22" + lat_dir + "%20" + to_string((int)lon_d) + "%C2%B0" + to_string((int)lon_m) + "%27" + to_string(lon_s) + "%22" + lon_dir;

    code = runExifTool("-GPS*", image_file, output_file);
    file = fopen(output_file.c_str(), "a");
    fprintf(file, "\n%s\n", maps_url.c_str());
    fclose(file);
    return 0;
}

string getImageFormat(const string filename)
{
    const string TEMP_FILE = "temp.txt";
    const int    BUFF_LEN  = 256;

    int code = runExifTool("-filetype", filename, TEMP_FILE);
    if (code != 0)
    {
        deleteFile(TEMP_FILE);
        return string("");
    }

    FILE* file = fopen(TEMP_FILE.c_str(), "r");
    char line[BUFF_LEN];
    char* ret = fgets(line, sizeof(line), file);
    fclose(file);
    deleteFile(TEMP_FILE);
    if (ret == NULL)
        return string("");

    int i = 0;
    char* token = strtok(line, " \n");
    while (token != NULL)
    {
        if (i == 3)
            return string(token);
        token = strtok(NULL, " \n");
        ++i;
    }
    return string("");
}

#endif // FILE_HPP