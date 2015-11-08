//
//  nodegen.h
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//  Copyright 2015. All rights reserved.
//
//

#ifndef __conSERN__nodegen__
#define __conSERN__nodegen__

#include <stdio.h>
#include "options.h"



typedef struct
{
    coordType *x;
    coordType *y;
    uint32_t *component;
} NodeList;

typedef struct
{
    double x;
    double y;
} VectorStruct; // maybe one day we will do 3d?


typedef struct
{
    uint32_t count;
    uint32_t allocated;
    VectorStruct *vertices;
} PolygonStruct;

typedef enum
{
    empty,
    partial,
    full
} BucketStatus;

typedef struct
{
    coordType    *x;
    coordType    *y;
    uint32_t count;
    uint32_t start;
    uint32_t i;
    uint32_t j;
    double  p;
    BucketStatus status;
    PolygonStruct *polygon;
} BucketStruct;


typedef enum
{
    rectangle,
    ellipse,
    polygon
} GeometryType;


typedef struct
{
    uint32_t Mx;
    uint32_t My;
    
    GeometryType type;
    VectorStruct origin;
    double xSize;
    double ySize;
    double totalArea;
    double bucketSize;
    PolygonStruct *polygon;
    
} GeometryStruct;

void * GenerateNodes(void *t);

void SetBucketSizes(const Options* options, const GeometryStruct* g, BucketStruct * buckets);

BucketStruct *GenerateBuckets(const Options *options, const GeometryStruct *g, NodeList *nodes);

GeometryStruct *geometryGenerate(const Options* options,
                                 GeometryType type,
                                 PolygonStruct *p);

void polygonAppend(const Options *options, PolygonStruct *p, const VectorStruct *v);

PolygonStruct *polygonNew(void);

void polygonFree(PolygonStruct *p);




#endif /* defined(__conSERN__nodegen__) */
