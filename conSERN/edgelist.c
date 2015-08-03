//
//  edgelist.c
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//
//
#include "FastSERN.h"
#include "edgelist.h"

uint32_t lock;
volatile uint64_t memccpyThreadCount;
pthread_mutex_t mutex_realloc;

uint64_t AllocateEdgeBuffer(EdgeList *l, uint32_t buffer_size, uint32_t enable_weights)
{
    l->weights_enabled = enable_weights;
    l->count = 0;
    l->allocated =  buffer_size;
    l->growth = 0; /* we don't use this member */
    
    
    l->from = malloc(sizeof(uint32_t) * l->allocated);
    l->to = malloc(sizeof(uint32_t) * l->allocated);
    
    
    if (l->weights_enabled)
        l->weight = malloc(sizeof(float) * l->allocated);
    
    return l->allocated;
}


void openEdgeList(EdgeList *l, Options* options)
{
    
    l->from = NULL;
    l->to = NULL;
    l->weight = NULL;
    l->component = NULL;
    l->weights_enabled = options->weights_enabled;
    l->count = 0;
    l->allocated = 0;
    l->growth = 0;
    
    pthread_mutex_init(&mutex_realloc, NULL);
    memccpyThreadCount = 0;
    lock =0;
}

void closeEdgeList(void)
{
    pthread_mutex_destroy(&mutex_realloc);
}

