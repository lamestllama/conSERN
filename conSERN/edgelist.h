//
//  edgelist.h
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//  Copyright 2015. All rights reserved.
//
//

#ifndef __conSERN__edgelist__
#define __conSERN__edgelist__

#include "options.h"
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>

extern uint64_t memccpyThreadCount;
extern pthread_mutex_t mutex_realloc;

typedef struct
{
    uint32_t *from;
    uint32_t *to;
    float    *weight;
    uint32_t weights_enabled;
    uint64_t count;
    uint64_t allocated;
    uint64_t growth;
} EdgeList;



uint64_t AllocateEdgeBuffer(Options *options, EdgeList *l, uint32_t buffer_size, uint32_t enable_weights);
void openEdgeList(EdgeList *l, Options *options);
void closeEdgeList(void);





// these functions are in the header because we want them to be inlined

static inline uint64_t CopyEdgeList(Options *options, EdgeList *l, EdgeList *b)
{
    uint64_t total;
    uint64_t count;
    
    
    // we need to decide if more memory is needed
    // then an how much inside a mutex the copy
    // could happen outside the mutex as long as it
    // doesn't rely on shared data.
    // The problem is what happens if some other thread
    // decides to do a reallocate while we are copying
    // the only thing our threads might mess up is stuff
    // pointed to by l
    
    pthread_mutex_lock(&mutex_realloc);
    
    count = l->count;
    total = l->count + b->count;
    
    /* if we need more memory */
    // TODO error handling for running out of memory or total too big
    if (total  >= l->allocated)
    {
        while (memccpyThreadCount != 0)
        {
            sched_yield();
        }
        
        /* if the growth amount is not enough extra */
        if ((l->allocated += l->growth) < total)
            l->allocated = total;
        
        /* realloc on a null pointer acts like malloc */
        l->from = options->realloc(l->from, sizeof(uint32_t) * l->allocated);
        l->to = options->realloc(l->to, sizeof(uint32_t) * l->allocated);
        
        /* only re-allocate space for the "distances" if required */
        if (l->weights_enabled)
            l->weight = options->realloc(l->weight,sizeof(float) * l->allocated);
        
        // if we needed memory but got none returned to us
        if ((l->allocated > 0) &&
            ((NULL == l->from) ||
             (NULL == l->to) ||
             (l->weights_enabled && (NULL == l->weight))))
        {
            
            options->errIdAndTxt(__FILE__, __LINE__,
                                 "Error unable to allocate memory");
        }
        
        /* TODO: we can do this better with some math */
        l->growth =  l->allocated / 10;
    }
    l->count = total;
    
    // atomic increment of count of write threads
    __sync_fetch_and_add(&memccpyThreadCount , 1 );
    
    // end critical section
    pthread_mutex_unlock(&mutex_realloc);
    
    
    memcpy(&(l->from[count]), b->from, b->count * sizeof(uint32_t));
    memcpy(&(l->to[count]), b->to, b->count * sizeof(uint32_t));
    
    if (l->weights_enabled)
        memcpy(&(l->weight[count]), b->weight,
               b->count * sizeof(float));
    
    b->count = 0;
    
    // atomic decrement of count of write threads.
    __sync_fetch_and_sub(&memccpyThreadCount , 1 );
    return l->allocated;
    
}


static inline uint64_t AddEdgeToBuffer(Options *options, EdgeList *l, EdgeList *b, uint32_t from, uint32_t to, double weight)
{
    
    if (b->count >= b->allocated)
    {
        CopyEdgeList(options, l, b);
    }
    
    
    b->from[b->count] = from;
    b->to[b->count] = to;
    
    if (b->weights_enabled) b->weight[b->count] = (float)weight;
    
    (b->count)++;
    
    return b->allocated;
}


#endif /* defined(__conSERN__edgelist__) */
