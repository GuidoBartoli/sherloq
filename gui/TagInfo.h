//------------------------------------------------------------------------------
// File:        TagInfo.h
//
// Description: Tag information object
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
#ifndef __TAGINFO_H__
#define __TAGINFO_H__

struct TagInfo
{
    TagInfo();
    virtual ~TagInfo();

    char    *group[3];  // family 0-2 group names
    char    *name;      // tag name
    char    *desc;      // tag description
    char    *id;        // tag ID
    char    *value;     // converted value
    int     valueLen;   // length of value in bytes (not including null terminator)
    char    *num;       // "numerical" value
    int     numLen;     // length of numerical value
    int     copyNum;    // copy number for this tag name
    TagInfo *next;      // next TagInfo in linked list
};

#endif // __TAGINFO_H__
