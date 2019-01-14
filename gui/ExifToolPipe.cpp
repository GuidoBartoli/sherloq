//------------------------------------------------------------------------------
// File:        ExifToolPipe.cpp
//
// Description: Piped output from exiftool application
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

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/wait.h>
#include "ExifToolPipe.h"

//------------------------------------------------------------------------------
// ExifToolPipe constructor
ExifToolPipe::ExifToolPipe()
            : mFile(-1), mBuff(NULL), mSize(0), mLen(0), mPos(0), mSearchPos(0),
              mBlockSize(0), mString(NULL), mStringLen(0), mPid(-1)
{
}

//------------------------------------------------------------------------------
// Destructor -- close the input file and delete the buffer
ExifToolPipe::~ExifToolPipe()
{
    if (mFile >= 0) close(mFile);
    Free();
}

//------------------------------------------------------------------------------
// Initialize the ExifTool pipe object
// - sets input file to non-blocking
void ExifToolPipe::Init(int fd, int pid, int initialSize)
{
    mFile = fd;
    mPid = pid;
    mBlockSize = initialSize;
    // set read pipe to non-blocking
    int flags = fcntl(fd, F_GETFL, 0) | O_NONBLOCK;
#ifdef O_BINARY
    // set to binary mode (Windows only)
    flags |= O_BINARY;
#endif
    fcntl(fd, F_SETFL, flags);
}

//------------------------------------------------------------------------------
// Read exiftool response for specified command number
// returns: command number on success, 0 if no response available, <0 on error
// - this routine returns immediately
int ExifToolPipe::Read()
{
    const int   kMinRemaining = 1024;   // enlarge buffer if less than this free

    Flush();    // remove previous response from buffer

    // keep reading until we get a complete response or there is no more to read
    for (;;) {

        // read to fill remaining response buffer, but leave room for null terminator
        int remaining = mSize - mLen - 1;

        // enlarge buffer if necessary
        // (could test for "remaining < 1", but what is the point in reading just 1 byte?)
        if (remaining < kMinRemaining) {
            int newSize = mSize + mLen + mBlockSize;
            char *pt = new char[newSize];
            if (!pt) return -3; // out of memory!
            if (mSize) memcpy(pt, mBuff, mSize);
            delete [] mBuff;
            mBuff = pt;
            mSize = newSize;
            remaining = newSize - mLen - 1;
        }

        // read output from exiftool process
        int bytesRead = (int)read(mFile, mBuff + mLen, remaining);
        if (bytesRead < 0) {
            if (errno != EAGAIN) return -2; // read error!
            bytesRead = 0;
        }
        mLen += bytesRead;
        if (mLen < 13) {
            if (!bytesRead) {
                // no response, so check to be sure our exiftool process is still running
                int status;
                if (mPid == -1) return -1;
                if (waitpid(mPid, &status, WNOHANG)) {
                    mPid = -1;
                    return -1;
                }
            }
            return 0;
        }
        // must null terminate response string
        mBuff[mLen] = '\0';
        //                           0123456789012
        // response should end with "{ready#####}\n" or "{ready#####}\r\n"
        // - continue searching from where we left off
        char *pt = mBuff + mSearchPos;
        char *end = mBuff + mLen;

        for (;;) {

            pt = (char *)memmem(pt, end-pt, "{ready", 6);
            if (!pt) {
                mSearchPos = mLen - 5;  // continue next search where we left off
                break;
            }
            if (end-pt>=13 && pt[11]=='}' &&
                // must end with newline, or CR+LF in Windows
                (pt[12]=='\n' || (pt[12]=='\r' && end-pt>=14 && pt[13]=='\n')))
            {
                // validate and extract command number
                int cmdNum = 0;
                for (int i=0; i<5; ++i) {
                    unsigned d = pt[i+6] -'0';
                    if (d > 9) {
                        cmdNum = 0;
                        break;
                    }
                    cmdNum = cmdNum * 10 + d;
                }
                if (cmdNum) {
                    *pt = '\0';         // NULL terminate original response
                    mStringLen = (int)(pt - mBuff);
                    // skip LF if this was a Windows CR+LF combo
                    if (pt[12] == '\r') ++pt;
                    pt += 13;           // step to start of next response
                    // update current position in response buffer
                    mPos = (int)(pt - mBuff);
                    mString = mBuff;    // set return string
                    return cmdNum;      // success!
                }
            }
            pt += 6;    // step to next possible search position
        }
        if (bytesRead != remaining) break;  // stop if we read everything
    }
    return 0;       // no complete response available
}

//------------------------------------------------------------------------------
// Free buffer memory
void ExifToolPipe::Free()
{
    delete [] mBuff;
    mBuff = mString = NULL;
    mLen = mSize = mPos = mSearchPos = mStringLen = 0;
}

//------------------------------------------------------------------------------
// Remove previous response from buffer
void ExifToolPipe::Flush()
{
    if (mPos) {
        if (mLen > mPos) {
            memmove(mBuff, mBuff+mPos, mLen-mPos);
            mLen -= mPos;
        } else {
            mLen = 0;
        }
        mPos = mSearchPos = 0;
    }
    mString = NULL;
    mStringLen = 0;
}

// end
