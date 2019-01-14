//------------------------------------------------------------------------------
// File:        ExifTool.h
//
// Description: C++ library interface to Perl exiftool application script
//
// License:     Copyright 2013-2016, Phil Harvey (phil at owl.phy.queensu.ca)
//
//              This is software, in whole or part, is free for use in
//              non-commercial applications, provided that this copyright notice
//              is retained.  A licensing fee may be required for use in a
//              commercial application.
//
// Created:     2013-11-23 - Phil Harvey
//------------------------------------------------------------------------------
#ifndef __EXIFTOOL_H__
#define __EXIFTOOL_H__

#include "ExifToolPipe.h"
#include "TagInfo.h"

#define NOW     0
#define NEVER   1e9

#define SUMMARY_DIRECTORIES_SCANNED     "directories scanned"
#define SUMMARY_DIRECTORIES_CREATED     "directories created"
#define SUMMARY_FILES_FAILED_CONDITION  "files failed condition"
#define SUMMARY_IMAGE_FILES_CREATED     "image files created"
#define SUMMARY_IMAGE_FILES_UPDATED     "image files updated"
#define SUMMARY_IMAGE_FILES_UNCHANGED   "image files unchanged"
#define SUMMARY_IMAGE_FILES_MOVED       "image files moved"
#define SUMMARY_IMAGE_FILES_COPIED      "image files copied"
#define SUMMARY_FILE_UPDATE_ERRORS      "files weren't updated due to errors"
#define SUMMARY_FILE_CREATE_ERRORS      "files weren't created due to errors"
#define SUMMARY_IMAGE_FILES_READ        "image files read"
#define SUMMARY_IMAGE_FILE_ERRORS       "files could not be read"
#define SUMMARY_OUTPUT_FILES_CREATED    "output files created"
#define SUMMARY_OUTPUT_FILES_APPENDED   "output files appended"
#define SUMMARY_HARD_LINKS_CREATED      "hard links created"
#define SUMMARY_HARD_LINK_ERRORS        "hard links could not be created"

class ExifTool
{
public:
            ExifTool(const char *exec=NULL, const char *arg1=NULL);
    virtual ~ExifTool();

    TagInfo *ImageInfo(const char *file, const char *opts=NULL, double timeout=NEVER);

    int     ExtractInfo(const char *file, const char *opts=NULL);
    TagInfo *GetInfo(int cmdNum=0, double timeout=NEVER);

    int     SetNewValue(const char *tag=NULL, const char *value=NULL, int len=-1);
    int     WriteInfo(const char *file, const char *opts=NULL, TagInfo *info=NULL);

    int     Command(const char *cmd=NULL);
    int     Complete(double timeout=NEVER);

    int     IsRunning();
    int     LastComplete()  { return mLastComplete; }
    int     LastCommand()   { return mCmdNum; } // (undocumented)
    void    SetLastComplete(int lastComplete) { mLastComplete = lastComplete; }

    char *  GetOutput()     { return mLastComplete > 0 ? mStdout.GetString() : NULL; }
    int     GetOutputLen()  { return mLastComplete > 0 ? mStdout.GetStringLen() : 0; }
    char *  GetError()      { return mLastComplete > 0 ? mStderr.GetString() : NULL; }
    int     GetErrorLen()   { return mLastComplete > 0 ? mStderr.GetStringLen() : 0; } // (undocumented)

    int     GetSummary(const char *msg);

    // flags to allow some ExifTool features to be disabled
    // (must be set before creating ExifTool object)
    static int  sNoSigPipe;     // set to disable SIGPIPE handler
    static int  sNoWatchdog;    // set to disable watchdog process

private:
    ExifToolPipe  mStdout;      // buffer for exiftool stdout read pipe
    ExifToolPipe  mStderr;      // buffer for exiftool stderr read pipe
    int           mTo;          // write pipe for exiftool stdin
    int           mPid;         // exiftool application process ID
    TagInfo     * mWriteInfo;   // tag information to write
    char        * mCmdQueue;    // queued command arguments (NULL if nothing queued)
    int           mCmdQueueLen; // length of data in command queue
    int           mCmdQueueSize;// size of command queue
    int           mWatchdog;    // watchdog process ID
    int           mLastComplete;// result of last Complete() call
    int           mCmdNum;      // last command number
};

#endif // __EXIFTOOL_H__
