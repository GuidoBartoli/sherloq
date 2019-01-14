#ifndef GENERAL_HPP
#define GENERAL_HPP

#include "utility.hpp"

using namespace std;

vector<string> automaticTagging(const string filename)
{
    string command = "python ../src/clarifai.py " + filename;
    string output = grabOutput(command);

    size_t start = output.find("{u'classes':") + 13;
    size_t end = output.find("]", start);
    string classes = output.substr(start, end - start);
    istringstream iss(classes);
    string token;
    vector<string> tags;
    while (getline(iss, token, ','))
    {
        string tag = token.substr(3, token.size() - 4);
        if (find(tags.begin(), tags.end(), tag) == tags.end())
            tags.push_back(tag);
    }

    return tags;
}

#endif // GENERAL_HPP