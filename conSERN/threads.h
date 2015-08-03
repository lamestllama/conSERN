//
//  threads.h
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//
//

#ifndef conSERN_threads_h
#define conSERN_threads_h

#include "options.h"
#include "edgelist.h"

typedef struct
{
    uint32_t     thread_id;
    uint32_t     thread_count;
    //uint32_t      N;
    //int          M;
    //double       s;
    //double       q;
    NodeList     *nodes;
    EdgeList     *edges;
    Options      *options;
    BucketStruct *buckets;
    GeometryStruct *geometry;
    double       *Q;
    uint32_t     BufferSize;
    
} ThreadDataStruct;

#endif
