//------------------------------------------------------------------------------
// File:        TagInfo.cpp
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

#include <stdlib.h>
#include "TagInfo.h"

//------------------------------------------------------------------------------
TagInfo::TagInfo()
       : name(NULL), desc(NULL), id(NULL), value(NULL), valueLen(0),
         num(NULL), numLen(0), copyNum(0), next(NULL)
{
    group[0] = group[1] = group[2] = NULL;
}

//------------------------------------------------------------------------------
// delete entire linked list of TagInfo objects
TagInfo::~TagInfo()
{
    // delete our  members
    delete [] group[0];
    delete [] group[1];
    delete [] group[2];
    delete [] name;
    delete [] desc;
    delete [] id;
    if (num != value) delete [] num;   // delete numerical value if unique
    delete [] value;

    // delete remaining elements of linked list
    while (next) {
        TagInfo *info = next;
        // remove next entry from the list, then delete it
        next = info->next;
        info->next = (TagInfo *)NULL;
        delete info;
    }
}

// end

