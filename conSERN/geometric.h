//
//  geometric.h
//  conSERN
//
//  Created by Eric Parsonage on 8/3/15.
//  Copyright 2015. All rights reserved.
//
//

#ifndef __conSERN__geometric__
#define __conSERN__geometric__

#include "uniform.h"
#include <math.h>

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif


static inline int64_t geom_rand(double p, uint32_t thread)
{
    // clamp the jumps to the maximum number of edges we can
    // have in a graph
    
    return (int64_t) min((double) INT64_MAX,
                          (1 + floor(log((GetUint(thread) + 1) *
                                         INTP_TO_DOUBLEP) / log(1 - p))));
    
   
}



/* lambda = -log2(1-p) */
static inline int64_t geom_rand2(double lambda, uint32_t thread)
{
    // clamp the jumps to the maximum number of edges we can
    // have in a graph addition of this will cause k to wrap negative
    // and then our main loop will exit the complexity in this method
    // is to avoid a branching instruction
    
    double results[2];
    
    results[0] = (-floor((log2(GetUint(thread) + 1) - 32) / lambda));
    results[1] = INT64_MAX;
    
    return (int64_t) results[results[0] >= results[1]];

   // and the alternative was this
   // return (int64_t) min((double) INT64_MAX,
   //                       (-floor((log2(GetUint(thread) + 1) - 32) / lambda)));
}

#endif /* defined(__conSERN__geometric__) */
