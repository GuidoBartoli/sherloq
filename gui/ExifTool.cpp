//------------------------------------------------------------------------------
// File:        ExifTool.cpp
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <fcntl.h>
#include <errno.h>
#include <ctype.h>
#include <sys/time.h>
#include <sys/wait.h>
#include "ExifTool.h"

const int   kOutBlockSize = 65536;  // size increment for exiftool stdout buffer
const int   kErrBlockSize = 4096;   // size increment for exiftool stderr buffer
const int   kCmdBlockSize = 8192;   // size increment for exiftool command buffer

const char *kDefaultExec = "exiftool";

static int  sBrokenPipe = -1;

int         ExifTool::sNoSigPipe = 0;
int         ExifTool::sNoWatchdog = 0;

//------------------------------------------------------------------------------
// SIGPIPE handler
static void sigPipeAction(int)
{
    sBrokenPipe = 1;
}

//------------------------------------------------------------------------------
// Unescape C-style escape sequences in a null-terminated string
// (as found in values of exiftool -php output)
// Returns: number of bytes in unescaped data (including null terminator)
// - on return string contains binary data which may have embedded zero bytes
static int unescape(char *str)
{
    char *src = strchr(str, '\\');
    // return string length without converting if string doesn't contain escapes
    if (!src) return((int)(strlen(str) + 1));
    char *dst = src;
    for (;;) {
        char ch = *(src++);
        if (ch == '\\') {
            // must un-escape this character
            ch = *(src++);
            switch (ch) {
                case 'x':
                    // decode 2-digit hex character
                    ch = 0;
                    for (int i=0; i<2; ++i) {
                        char nibble = *(src++);
                        if (nibble >= '0' && nibble <= '9') {
                            nibble -= '0';
                        } else if (nibble >= 'A' && nibble <= 'F') {
                            nibble -= 'A' - 10;
                        } else if (nibble >= 'a' && nibble <= 'f') {
                            nibble -= 'a' - 10;
                        } else {
                            ch = 0; // (shouldn't happen)
                            break;
                        }
                        ch = (ch << 4) + nibble;
                    }
                    break;
                case 't':
                    ch = '\t';
                    break;
                case 'n':
                    ch = '\n';
                    break;
                case 'r':
                    ch = '\r';
                    break;
                case '\0':  // (shouldn't happen, but just to be safe)
                    *(dst++) = ch;
                    return((int)(dst - str));
                default:
                    // pass any other character straight through
                    break;
            }
            *(dst++) = ch;
        } else {
            *(dst++) = ch;
            if (!ch) break;
        }
    }
    return((int)(dst - str));
}

//------------------------------------------------------------------------------
// get current hi-resolution time
static double getTime()
{
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv,&tz);
    return(tv.tv_sec + 1e-6 * tv.tv_usec);
}

//==============================================================================
// ExifTool object constructor
//------------------------------------------------------------------------------
// Inputs: exec - path name to executable file (ie. "perl" or "exiftool")
//         arg1 - optional first argument (ie. "exiftool" if exec="perl")
ExifTool::ExifTool(const char *exec, const char *arg1)
        : mWriteInfo(NULL), mCmdQueue(NULL), mCmdQueueLen(0), mCmdQueueSize(0),
          mWatchdog(-1), mLastComplete(0), mCmdNum(0)
{
    int to[2], from[2], err[2];
    const char *args[7];
    args[0] = NULL;
    args[1] = kDefaultExec;
    args[2] = "-stay_open";
    args[3] = "true";
    args[4] = "-@";
    args[5] = "-";
    args[6] = NULL;

    int firstArg = 1;
    if (arg1) {
        args[1] = arg1;
        --firstArg;
    }
    args[firstArg] = exec ? exec : kDefaultExec;

    // set up handler for SIGPIPE if not done already
    if (sBrokenPipe == -1 && !sNoSigPipe) {
        struct sigaction act1;
        sigemptyset(&act1.sa_mask);
        act1.sa_flags = SA_SIGINFO;
        act1.sa_handler =  sigPipeAction;
        if (sigaction(SIGPIPE, &act1, (struct sigaction *)NULL) < 0) {
            sBrokenPipe = -2;   // error!
        } else {
            sBrokenPipe = 0;
        }
    }

    // create our pipes
    pipe(to);
    pipe(from);
    pipe(err);

    // fork and exec exiftool
    mPid = fork();

    if (mPid == 0) {
        // this is our exiftool thread
        close(to[1]);
        close(from[0]);
        close(err[0]);
        dup2(to[0], STDIN_FILENO);
        dup2(from[1], STDOUT_FILENO);
        dup2(err[1], STDERR_FILENO);
        close(to[0]);
        close(from[1]);
        close(err[1]);
        execvp(args[firstArg], (char * const *)args + firstArg);
        // (if execvp succeeds, it will never return)
        exit(0);
    } else {
        // this is our working thread
        close(to[0]);
        close(from[1]);
        close(err[1]);
        // allocate memory for exiftool response and error messages
        mStdout.Init(from[0], mPid, kOutBlockSize);
        mStderr.Init(err[0], mPid, kErrBlockSize);
        mTo = to[1];
        // set output exiftool argument pipe to non-blocking
        int flags = fcntl(mTo, F_GETFL, 0);
        fcntl(mTo, F_SETFL, flags | O_NONBLOCK);
        // create a watchdog process to monitor the main program thread and
        // terminate the exiftool process if necessary when the program exits
        // (otherwise exiftool processes may be left running if an ExifTool
        // object wasn't properly deleted)
        if (mPid != -1 && !sNoWatchdog) {
            int pid = getpid(); // pid of this process
            mWatchdog = fork();
            if (!mWatchdog) {
                // monitor the main program and clean up if it dies
                // (under normal conditions this watchdog process
                // is killed in the destructor)
                for (;;) {
                    sleep(1);
                    // if the main thread dies, our parent PID will change
                    if (getppid() == pid) continue;
                    // our parent has died, so send a SIGINT to exiftool
                    kill(mPid, SIGINT);
                    exit(0);    // exit the watchdog process
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
// Terminate exiftool application process and close files
ExifTool::~ExifTool()
{
    delete mWriteInfo;
    delete [] mCmdQueue;
    Command("-stay_open\nfalse\n");
    close(mTo);
    // wait for the exiftool process to terminate
    while (IsRunning()) {
        // keep reading exiftool output so the process won't be blocked,
        // and wait until all commands complete unless we get an error
        if (Complete() < 0) {
            kill(mPid, SIGINT);     // pull the plug
            break;
        }
        usleep(100);    // relax
    }
    if (mWatchdog > 0) kill(mWatchdog, SIGINT); // kill our watchdog process
}

//------------------------------------------------------------------------------
// Extract metadata from specified image
// Inputs:  file - source file name
//          opts - string of exiftool options, separated by newlines
//          timeout - maximum wait time (floating point seconds)
// Returns: linked list of tag information structures for extracted information
// - valid options: "-b\n" - extract binary data (other options may mess up TagInfo parsing)
// - waits for exiftool command to complete, then parses returned messages
//   to generate the TagInfo list
// - caller is responsible for deleting returned TagInfo ("delete info;")
TagInfo * ExifTool::ImageInfo(const char *file, const char *opts, double timeout)
{
    int cmdNum = ExtractInfo(file, opts);

    // error unless command number is > 0
    if (cmdNum <= 0) return (TagInfo *)NULL;

    return GetInfo(cmdNum, timeout);
}

//------------------------------------------------------------------------------
// Send command to extract information from one or more files
// Inputs: file - source file name(s)
//         opts - exiftool options, separated by newlines
// Returns: command number (>0), or error (<0)
int ExifTool::ExtractInfo(const char *file, const char *opts)
{
    if (!file) return -5;

    // prepare command arguments for exiftool
    int flen = (int)strlen(file);
    int olen = opts ? (int)strlen(opts) : 0;

    char *buff = new char[flen + olen + 64];
    if (!buff) return -3;       // out of memory!
    memcpy(buff, file, flen);

    // extract information using exiftool -php option
    // (also add "-l -G:0:1:2:4 -D" to get more details for TagInfo structure)
    strcpy(buff + flen, "\n-php\n-l\n-G:0:1:2:4\n-D\n-sep\n, \n");

    // add extra options specified by the caller
    if (opts) {
        strcat(buff, opts);
        strcat(buff, "\n"); // (to be safe; blank lines are ignored)
    }

    // send the command to exiftool
    int cmdNum = Command(buff);
    delete [] buff;

    return cmdNum;
}

//------------------------------------------------------------------------------
// Read exiftool output and convert to list of TagInfo structures
// Inputs:  cmdNum - command number (0 to process next output in series,
//                      or -1 to process previously completed command)
//          timeout - maximum wait time (floating point seconds)
// Returns: linked list of tag information structures for extracted information
// - waits up to timeout time for exiftool command to complete
// - caller is responsible for deleting returned TagInfo ("delete info;")
// - after return, GetError() may be called to get exiftool errors
TagInfo *ExifTool::GetInfo(int cmdNum, double timeout)
{
    TagInfo *info = (TagInfo *)NULL;
    TagInfo *next = (TagInfo *)NULL;
    TagInfo *infoList = (TagInfo *)NULL;

    // wait for specified command to complete
    if (cmdNum >= 0) {
        for (;;) {
            int n = Complete(timeout);
            if (n <= 0) return info;
            if (n == cmdNum || !cmdNum) break;
        }
    } else if (mLastComplete <= 0) {
        return info;
    }

    // parse the response string, line by line
    char *pt = mStdout.GetString();
    if (!pt) return info;

    int mode = 0;   // 0=looking for tag name, 1=tag properties

    for (;;) {
        // find the end of this line
        char *p0 = pt;
        pt = strchr(pt, '\n');
        if (!pt) break;
        *pt = '\0'; // null terminate this line
        ++pt;       // continue at next line
        // scan for opening quote
        p0 = strchr(p0, '"');
        if (!p0) {
            // finish with most recent TagInfo structure
            if (info) {
                // (name and value are guaranteed to exist)
                if (!info->value) {
                    info->value = new char[1];
                    if (!info->value) break;
                    info->value[0] = '\0';
                }
                if (!info->num) {
                    info->num = info->value;
                    info->numLen = info->valueLen;
                }
            }
            mode = 0;               // look for next tag name
            continue;
        }
        char *p1 = ++p0;
        if (!mode) {    // looking for new tag
            if (next) delete next;  // delete old unused structure if necessary
            // create new TagInfo structure for this tag
            next = new TagInfo;
            if (!next) break;
            // extract tag/group names
            int g = 0;
            for (;;) {
                char ch = *p1;
                if (ch == '"' || ch == ':') {
                    int n = (int)(p1 - p0);
                    char *str = new char[n + 1];
                    if (!str) break;
                    memcpy(str, p0, n);
                    str[n] = '\0';
                    if (ch == '"') {
                        next->name = str;   // save tag name
                        break;
                    }
                    if (g > 2) {
                        // get copy number
                        if (!memcmp(str, "Copy", 4)) {
                            next->copyNum = atoi(str+4);
                            delete [] str; // done with this string
                        }
                    } else {
                        next->group[g] = str;   // save group name
                    }
                    ++g;
                    p0 = p1 + 1;
                }
                ++p1;
            }
            if (!next->name) continue;
            // file name given by line like:  "SourceFile" => "images/a.jpg",
            if (!strcmp(next->name,"SourceFile")) {
                char *p2 = pt - 2;
                if (*p2 == '\r') --p2; // skip Windows CR
                if (*p2 == ',') --p2;
                if (*p2 != '"') continue;
                int n = (int)(p2 - p1 - 6);
                if (n < 0) continue;
                char *str = new char[n+1];
                if (!str) break;
                memcpy(str, p1+6, n);
                str[n] = '\0';
                next->value = str;
                next->valueLen = n;
            } else {
                mode = 1;   // read tag properties next
            }
            // add to linked list of information
            if (info) {
                info->next = next;
            } else {
                infoList = next;
            }
            info = next;
            next = NULL;
        } else {
            // isolate the property name
            p1 = strchr(p0, '"');
            if (!p1) break;         // (shouldn't happen)
            *p1 = '\0';             // null terminate property name
            p1 += 5;                // step to start of value
            if (p1 >= pt) break;    // (shouldn't happen);
            if (*p1 == '"') ++p1;   // skip quote if it exists
            char *p2 = pt - 1;
            if (p2[-1] == '\r') --p2;// skip Windows CR
            if (p2[-1] == ',') --p2;// skip trailing comma
            if (p2[-1] == '"') --p2;// skip trailing quote
            if (p2 < p1) break;     // (shouldn't happen)
            *p2 = '\0';             // null terminate property value
            int n = unescape(p1);   // unescape characters in property value
            char **dst;
            if (!strcmp(p0, "desc")) {
                dst = &info->desc;
            } else if (!strcmp(p0, "id")) {
                dst = &info->id;
            } else if (!strcmp(p0, "num")) {
                dst = &info->num;
                info->numLen = n - 1;   // save length too (could be binary data)
            } else if (!strcmp(p0, "val")) {
                dst = &info->value;
                info->valueLen = n - 1; // save length too (could be binary data)
            } else {
                continue;   // (shouldn't happen)
            }
            *dst = new char[n];
            if (!*dst) break;
            memcpy(*dst, p1, n);    // copy property value
        }
    }
    if (next) delete next;
    return infoList;
}

//------------------------------------------------------------------------------
// Set the new value for a tag
// Inputs:  tag = tag name (may contain leading group names and trailing '#')
//          value = tag value data
//          len = length of value in bytes (defaults to strlen(value))
// Returns: number of tags set, or <0 on memory error
// - must call WriteInfo() at some point after this to actually write the new values
// - call with tag=NULL to reset previous new values
// - call with value=NULL to delete tag
int ExifTool::SetNewValue(const char *tag, const char *value, int len)
{
    int numSet = 0;
    if (tag) {
        TagInfo *info = new TagInfo;
        if (!info) return -3;
        info->name = new char[strlen(tag) + 1];
        if (!info->name) { delete info; return -3; }
        strcpy(info->name, tag);
        if (value) {
            if (len < 0) {
                if (value) len = (int)strlen(value);
                else len = 0;
            }
            if (len) {
                info->value = new char[len+1];
                if (!info->value) { delete info; return -3; }
                memcpy(info->value, value, len);
                // add null terminator (but note that we don't use it)
                info->value[len] = '\0';
                info->valueLen = len;
            }
        }
        // place at the end of the linked list
        TagInfo **pt = &mWriteInfo;
        while (*pt) {
            ++numSet;
            pt = &((*pt)->next);
        }
        *pt = info;
        ++numSet;
    } else {
        delete mWriteInfo;
        mWriteInfo = NULL;
    }
    return numSet;
}

//------------------------------------------------------------------------------
// Write metadata to specified image
// Inputs:  file - one or more directory and/or file names, separated by newlines
//          opts - extra exiftool options, separated by newlines
//          info - pointer to linked list of tags to write (overrides SetNewValue() calls)
// Returns: >0 command number, <0=error
// - each option argument must be terminated with a newline
// - file may be one or more file names separated by newlines
// - must call Complete() after this to check the return messages
// - ignores "SourceFile" entries in input TagInfo list
int ExifTool::WriteInfo(const char *file, const char *opts, TagInfo *info)
{
    if (!file) return -5;

    const int kBlockSize = 65536;
    char *buff = new char[kBlockSize];
    if (!buff) return -3;
    // get length of all options (plus 12 extra characters for "-ex\n-sep\n, \n")
    int olen = (int)((opts ? strlen(opts)+1 : 0) + 12);
    int size = kBlockSize;
    strcpy(buff, file);
    int pos = (int)strlen(file);
    buff[pos++] = '\n';
    int escaped = 0;

    // use internal new value list if not specified
    if (!info) info = mWriteInfo;

    // prepare assignment arguments for exiftool, looping through all tags to write
    while (info) {
        if (!info->name || strlen(info->name) > 100 || !strcmp(info->name, "SourceFile")) {
            info = info->next;
            continue;
        }
        // construct the tag name
        char tag[1024];
        tag[0] = '\0';
        for (int i=0; i<3; ++i) {
            if (info->group[i] && strlen(info->group[i]) < 100) {
                char *pt = strchr(tag, '\0');
                *(pt++) = '0' + i;          // leading family number
                strcpy(pt, info->group[i]); // group name
                strcat(tag,":");            // colon separator
            }
        }
        strcat(tag, info->name);
        // which value are we writing (converted or numerical?)
        char *val = info->value;
        int valLen = info->valueLen;
        if (!val) {
            val = info->num;
            valLen = info->numLen;
            if (val) strcat(tag, "#");  // write numerical value
        }
        int tagLen = (int)strlen(tag);
        int origLen = valLen;
        // count the number of characters in the value that need escaping
        if (val) {
            char *pt = val;
            char *end = pt + origLen;
            int count = 0;
            while (pt < end) {
                char ch = *(pt++);
                if (ch==10 || ch=='&' || ch=='\0') ++count;
            }
            valLen += count * 4;   // 4 extra bytes for each escaped character
        }
        // prepare exiftool argument (format is: "-TAG=VALUE\n")
        int n = tagLen + valLen + 3;
        // expand buffer if necessary
        if (pos + n + olen > size) {
            size += n + kBlockSize;
            char *buf2 = new char[size];
            if (!buf2) { delete [] buff; return -3; }
            memcpy(buf2, buff, pos);
            delete [] buff;
            buff = buf2;
        }
        buff[pos++] = '-';
        memcpy(buff+pos, tag, tagLen);
        pos += tagLen;
        buff[pos++] = '=';
        if (valLen == origLen) {
            if (val) {
                memcpy(buff+pos, val, valLen);
                pos += valLen;
            }
        } else {
            // escape newlines and '&' characters in the value
            char *pt = val;
            char *end = pt + origLen;
            char *dst = buff + pos;
            while (pt < end) {
                char ch = *(pt++);
                if (ch == 10) {
                    memcpy(dst, "&#10;", 5);
                    dst += 5;
                } else if (ch == 0) {
                    memcpy(dst, "&#00;", 5);
                    dst += 5;
                } else if (ch == '&') {
                    memcpy(dst, "&amp;", 5);
                    dst += 5;
                } else {
                    *(dst++) = ch;
                }
            }
            pos = (int)(dst - buff);
            escaped = 1;
        }
        buff[pos++] = '\n';
        info = info->next;
    }
    // get exiftool to unescape our newlines if necessary
    if (escaped) { strcpy(buff+pos, "-ex\n");  pos += 4; }
    // split concatenated lists back into individual items
    strcpy(buff+pos, "-sep\n, \n");  pos += 8;
    // add user-defined options last
    if (opts) {
        strcpy(buff+pos, opts);
        pos += strlen(opts);
        buff[pos++] = '\n'; // (just in case)
        buff[pos] = '\0';
    }
    int cmdNum = Command(buff);
    delete [] buff;
    return cmdNum;
}

//------------------------------------------------------------------------------
// Send command to exiftool
// Inputs:  cmd - exiftool command arguments (separated by newlines)
// Returns: command number (1-99999), error (<0),
//          or 0 if cmd is NULL and the queue is empty
// - set cmd to NULL to continue writing previously queued commands
// - this routine always returns immediately, and queues command if it couldn't send
int ExifTool::Command(const char *cmd)
{
    int n;
    // check to make sure our exiftool process is still running
    if (!IsRunning()) return -1;
    // must first try to send previously queued command
    if (mCmdQueue) {
        n = (int)write(mTo, mCmdQueue, mCmdQueueLen);
        if (n < 0) {
            if (errno != EAGAIN) return -2; // write error!
            n = 0;
        }
        if (n == mCmdQueueLen) {
            delete [] mCmdQueue;
            mCmdQueue = NULL;
            mCmdQueueSize = 0;
        } else if (n != 0) {
            memmove(mCmdQueue, mCmdQueue+n, mCmdQueueLen-n);
        }
        mCmdQueueLen -= n;
    }
    if (cmd) {
        // increment to next command number
        int cmdNum = mCmdNum + 1;
        if (cmdNum > 99999) cmdNum = 1;

        // compose full command string (cmd2)
        char buf2[64];
        int len = (int)strlen(cmd);
        int len2 = sprintf(buf2, "\n-echo4\n{ready%.5d}\n-execute%.5d\n", cmdNum, cmdNum);
        char *cmd2 = new char[len + len2];
        if (!cmd2) return -3;        // out of memory!
        memcpy(cmd2, cmd, len);
        memcpy(cmd2+len, buf2, len2);
        len2 += len;    // len2 is now the length of the complete command

        // add command to the queue if not empty
        if (mCmdQueue) {
            if (mCmdQueueLen + len2 > mCmdQueueSize) {
                // enlarge queue and add new command
                int newSize = mCmdQueueLen + len2 + kCmdBlockSize;
                char *queue = new char[newSize];
                if (!queue) return -3;          // out of memory!
                memcpy(queue, mCmdQueue, mCmdQueueLen);
                delete [] mCmdQueue;
                mCmdQueue = queue;
                mCmdQueueSize = newSize;
            }
            // copy this command into the queue
            memcpy(mCmdQueue+mCmdQueueLen, cmd2, len2);
            mCmdQueueLen += len2;
            delete [] cmd2;     // free memory for this command
        } else {
            // write the command
            n = (int)write(mTo, cmd2, len2);
            if (n < 0) {
                if (errno != EAGAIN) return -2; // write error!
                n = 0;
            }
            if (n == len2) {
                delete [] cmd2; // success! delete the buffered command
            } else {
                // don't bother allocating any new memory,
                // just use our cmd2 string as the new queue
                if (n) memmove(cmd2, cmd2+n, len2-n);
                mCmdQueue = cmd2;
                mCmdQueueLen = len2 - n;
                mCmdQueueSize = len2;
            }
        }
        mCmdNum = cmdNum;
    } else if (!mCmdQueue) {
        return 0;   // cmd is NULL, and queue is empty
    }
    return mCmdNum;
}

//------------------------------------------------------------------------------
// Wait for a command to complete (up to specified timeout)
// Inputs:  timeout - maximum wait time (floating point seconds)
// Returns: command number on success, 0 on timeout, or <0 on error
int ExifTool::Complete(double timeout)
{
    if (mCmdQueue) Command();       // try to send queued commands (if any)
    double doneTime = getTime() + timeout;
    int cmdNum;
    for (;;) {
        cmdNum = mStdout.Read();
        if (cmdNum) break;
        if (getTime() >= doneTime) break;
        if (mCmdQueue) Command();   // keep sending queued commands
        usleep(1000);               // chill and have a beer
    }
    if (cmdNum > 0) {
        // get errors from the same command (we know they must be coming,
        // so loop as quickly as possible to read them, but impose a
        // 1-second timeout just in case)
        doneTime = getTime() + 1;
        for (;;) {
            int n = mStderr.Read();
            if (n == cmdNum) break;
            if (n < 0) { cmdNum = n; break; }
            if (getTime() >= doneTime) { cmdNum = -4; break; }
        }
    }
    return mLastComplete = cmdNum;
}

//------------------------------------------------------------------------------
// Get specified summary message
// Inputs: msg - message string in summary output
// Returns: corresponding number from summary statistics, or -1 if the
//          specified message wasn't found
int ExifTool::GetSummary(const char *msg)
{
    for (int out=0; out<2; ++out) {
        // check stderr first because it isn't likely to be too long
        char *str = out ? GetOutput() : GetError();
        if (!str) continue;
        char *pt = strstr(str, msg);
        if (!pt || pt - str < 2 || pt[-1] != ' ' || !isdigit(pt[-2])) continue;
        char ch = pt[strlen(msg)];
        if (ch != '\n' && ch != '\r') continue; // message must end with a newline
        pt -= 2;
        while (pt > str && isdigit(pt[-1])) --pt;
        return atoi(pt);
    }
    return -1;  // message not found
}

//------------------------------------------------------------------------------
// Check to see if exiftool process is still running
int ExifTool::IsRunning()
{
    int status;
    if (mPid == -1) return 0;
    if (waitpid(mPid, &status, WNOHANG)) {
        // no more child process
        mPid = -1;
        return 0;
    }
    return 1;   // yes!
}

// end

