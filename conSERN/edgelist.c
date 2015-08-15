//
//  edgelist.c
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//  Copyright 2015. All rights reserved.
//
//
#include "FastSERN.h"
#include "edgelist.h"

uint32_t lock;
uint64_t memccpyThreadCount;
pthread_mutex_t mutex_realloc;

uint64_t AllocateEdgeBuffer(Options *options, EdgeList *l, uint32_t buffer_size, uint32_t enable_weights)
{
    l->weights_enabled = enable_weights;
    l->count = 0;
    l->allocated =  buffer_size;
    l->growth = 0; /* we don't use this member */
    
    
    l->from = malloc(sizeof(uint32_t) * l->allocated);
    if (l->from == NULL)
        options->errIdAndTxt(__FILE__, __LINE__,
                             "Error unable to allocate memory");
    
    l->to = malloc(sizeof(uint32_t) * l->allocated);
    if (l->to == NULL)
        options->errIdAndTxt(__FILE__, __LINE__,
                             "Error unable to allocate memory");
    
    
    if (l->weights_enabled)
    {
        l->weight = malloc(sizeof(float) * l->allocated);
        if (l->weight == NULL)
            options->errIdAndTxt(__FILE__, __LINE__,
                                 "Error unable to allocate memory");
    }
    
    return l->allocated;
}


void openEdgeList(EdgeList *l, Options* options)
{
    memset(l, 0, sizeof(EdgeList));
   
    l->weights_enabled = options->weights_enabled;
   
    pthread_mutex_init(&mutex_realloc, NULL);
    memccpyThreadCount = 0;
    lock =0;
}

void closeEdgeList(void)
{
    pthread_mutex_destroy(&mutex_realloc);
}

